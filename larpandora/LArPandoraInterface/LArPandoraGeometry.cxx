/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.cxx
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#include "cetlib/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "Pandora/PandoraInternal.h" // For pandora::TypeToString
#include "Xml/tinyxml.h"

#include <iomanip>
#include <sstream>

namespace lar_pandora
{

void LArPandoraGeometry::LoadDetectorGaps(const Settings &settings, LArDetectorGapList &listOfGaps)
{
    // Detector gaps can only be loaded once - throw an exception if the output lists are already filled
    if (!listOfGaps.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadDetectorGaps --- the list of gaps already exists ";

    // Only populate detector gaps when using global drift volume
    if (!settings.m_globalDriftVolume)
        return;

    // Loop over drift volumes and write out the dead regions at their boundaries
    LArDriftVolumeList driftVolumeList;
    LArPandoraGeometry::LoadGeometry(driftVolumeList);

    for (LArDriftVolumeList::const_iterator iter1 = driftVolumeList.begin(), iterEnd1 = driftVolumeList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArDriftVolume &driftVolume1 = *iter1;

        for (LArDriftVolumeList::const_iterator iter2 = iter1, iterEnd2 = driftVolumeList.end(); iter2 != iterEnd2; ++iter2)
        {
            const LArDriftVolume &driftVolume2 = *iter2;

            if (driftVolume1.GetVolumeID() == driftVolume2.GetVolumeID())
                continue;

            const double maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?
            const double deltaZ(std::fabs(driftVolume1.GetCenterZ() - driftVolume2.GetCenterZ()));
            const double deltaY(std::fabs(driftVolume1.GetCenterY() - driftVolume2.GetCenterY()));
            const double deltaX(std::fabs(driftVolume1.GetCenterX() - driftVolume2.GetCenterX()));
            const double widthX(0.5f * (driftVolume1.GetWidthX() + driftVolume2.GetWidthX()));
            const double gapX(deltaX - widthX);

            if (gapX < 0.0 || gapX > maxDisplacement || deltaY > maxDisplacement || deltaZ > maxDisplacement)
                continue;

            const double X1((driftVolume1.GetCenterX() < driftVolume2.GetCenterX()) ? (driftVolume1.GetCenterX() + 0.5 * driftVolume1.GetWidthX()) :
                (driftVolume2.GetCenterX() + 0.5 * driftVolume2.GetWidthX()));
            const double X2((driftVolume1.GetCenterX() > driftVolume2.GetCenterX()) ? (driftVolume1.GetCenterX() - 0.5 * driftVolume1.GetWidthX()) :
                (driftVolume2.GetCenterX() - 0.5 * driftVolume2.GetWidthX()));
            const double Y1(std::min((driftVolume1.GetCenterY() - 0.5 * driftVolume1.GetWidthY()),
                (driftVolume2.GetCenterY() - 0.5 * driftVolume2.GetWidthY())));
            const double Y2(std::max((driftVolume1.GetCenterY() + 0.5 * driftVolume1.GetWidthY()),
                (driftVolume2.GetCenterY() + 0.5 * driftVolume2.GetWidthY())));
            const double Z1(std::min((driftVolume1.GetCenterZ() - 0.5 * driftVolume1.GetWidthZ()),
                (driftVolume2.GetCenterZ() - 0.5 * driftVolume2.GetWidthZ())));
            const double Z2(std::max((driftVolume1.GetCenterZ() + 0.5 * driftVolume1.GetWidthZ()),
                (driftVolume2.GetCenterZ() + 0.5 * driftVolume2.GetWidthZ())));

            listOfGaps.push_back(LArDetectorGap(X1, Y1, Z1, X2, Y2, Z2));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::LoadGeometry(const Settings &settings, LArDriftVolumeList &outputVolumeList, LArDriftVolumeMap &outputVolumeMap)
{
    // Geometry can only be loaded once - throw an exception if the output lists are already filled
    if (!outputVolumeList.empty() || !outputVolumeMap.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- the list of drift volumes already exists ";

    // Global drift volume - use a global coordinate system and combine all drift volumes into a single global volume
    if (settings.m_globalDriftVolume)
    {
        LArDriftVolumeList inputVolumeList;
        LArPandoraGeometry::LoadGeometry(inputVolumeList);
        LArPandoraGeometry::LoadGlobalParentGeometry(inputVolumeList, outputVolumeList);
    }

    // Global coordinates - use a global coordinate system but keep drift volumes separate
    else if (settings.m_globalCoordinates)
    {
        LArDriftVolumeList inputVolumeList;
        LArPandoraGeometry::LoadGeometry(inputVolumeList);
        LArPandoraGeometry::LoadGlobalDaughterGeometry(inputVolumeList, outputVolumeList);
    }

    // Separate drift volumes - keep drift volumes separate with their separate coordinate systems
    else
    {
        LArPandoraGeometry::LoadGeometry(outputVolumeList);
    }

    // Create mapping between tpc/cstat labels and drift volumes
    for (const LArDriftVolume &driftVolume : outputVolumeList)
    {
        for (const LArTpcVolume &tpcVolume : driftVolume.GetTpcVolumeList())
        {
            (void) outputVolumeMap.insert(LArDriftVolumeMap::value_type(LArPandoraGeometry::GetTpcID(tpcVolume.GetCryostat(), tpcVolume.GetTpc()), driftVolume));
        }
    }

    // Print drift volumes for debugging
    if (settings.m_printGeometry)
        LArPandoraGeometry::PrintGeometry(settings, outputVolumeList);

    if (!settings.m_geometryXmlFileName.empty())
        LArPandoraGeometry::WriteGeometry(settings.m_geometryXmlFileName, outputVolumeList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::PrintGeometry(const Settings &settings, const LArDriftVolumeList &driftVolumeList)
{
    // Useful for debugging (or reverse-engineering) a detector geometry
    std::cout << " *** LArPandoraGeometry::PrintGeometry(...) *** " << std::endl;

    std::cout << " *** Geometry Settings *** " << std::endl;
    std::cout << "  globalCoordinates = " << settings.m_globalCoordinates << std::endl;
    std::cout << "  globalDriftVolume = " << settings.m_globalDriftVolume << std::endl;

    for (const LArDriftVolume &driftVolume : driftVolumeList)
    {
        std::cout << " *** Drift Volume *** " << std::endl;
        std::cout << "  volumeID = " << driftVolume.GetVolumeID() << std::endl;
        std::cout << "  isPositiveDrift = " << driftVolume.IsPositiveDrift() << std::endl;
        std::cout << "  wirePitchU = " << driftVolume.GetWirePitchU() << std::endl;
        std::cout << "  wirePitchV = " << driftVolume.GetWirePitchV() << std::endl;
        std::cout << "  wirePitchW = " << driftVolume.GetWirePitchW() << std::endl;
        std::cout << "  wireAngleU = " << driftVolume.GetWireAngleU() << std::endl;
        std::cout << "  wireAngleV = " << driftVolume.GetWireAngleV() << std::endl;
        std::cout << "  centerX = " << driftVolume.GetCenterX() << std::endl;
        std::cout << "  centerY = " << driftVolume.GetCenterY() << std::endl;
        std::cout << "  centerZ = " << driftVolume.GetCenterZ() << std::endl;
        std::cout << "  widthX = " << driftVolume.GetWidthX() << std::endl;
        std::cout << "  widthY = " << driftVolume.GetWidthY() << std::endl;
        std::cout << "  widthZ = " << driftVolume.GetWidthZ() << std::endl;
        std::cout << "  sigmaUVZ = " << driftVolume.GetSigmaUVZ() << std::endl;

        std::cout << "  TPC LIST [cstat][tpc]: " << std::endl;

        for (const LArTpcVolume &tpcVolume : driftVolume.GetTpcVolumeList())
            std::cout << "   [" << tpcVolume.GetCryostat() << "][" << tpcVolume.GetTpc() << "]" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::WriteGeometry(const std::string &xmlFileName, const LArDriftVolumeList &driftVolumeList)
{
    pandora::TiXmlDocument xmlDocument;

    for (const LArDriftVolume &driftVolume : driftVolumeList)
    {
        pandora::TiXmlElement *const pVolumeElement = new pandora::TiXmlElement("LArDriftVolume");

        pandora::TiXmlElement *const pVolumeIdElement = new pandora::TiXmlElement("VolumeID");
        pVolumeIdElement->LinkEndChild(new pandora::TiXmlText(pandora::TypeToString(driftVolume.GetVolumeID())));
        pVolumeElement->LinkEndChild(pVolumeIdElement);

        pandora::TiXmlElement *const pDriftDirectionElement = new pandora::TiXmlElement("IsPositiveDrift");
        pDriftDirectionElement->LinkEndChild(new pandora::TiXmlText(pandora::TypeToString(driftVolume.IsPositiveDrift())));
        pVolumeElement->LinkEndChild(pDriftDirectionElement);

        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WirePitchU", driftVolume.GetWirePitchU());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WirePitchV", driftVolume.GetWirePitchV());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WirePitchW", driftVolume.GetWirePitchW());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WireAngleU", driftVolume.GetWireAngleU());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WireAngleV", driftVolume.GetWireAngleV());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "CenterX", driftVolume.GetCenterX());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "CenterY", driftVolume.GetCenterY());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "CenterZ", driftVolume.GetCenterZ());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WidthX", driftVolume.GetWidthX());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WidthY", driftVolume.GetWidthY());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "WidthZ", driftVolume.GetWidthZ());
        LArPandoraGeometry::WritePrecisionElement(pVolumeElement, "SigmaUVZ", driftVolume.GetSigmaUVZ());

        xmlDocument.LinkEndChild(pVolumeElement);
    }

    xmlDocument.SaveFile(xmlFileName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPandoraGeometry::GetVolumeID(const LArDriftVolumeMap &driftVolumeMap, const unsigned int cstat, const unsigned int tpc)
{
    if (driftVolumeMap.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::GetVolumeID --- detector geometry map is empty";

    LArDriftVolumeMap::const_iterator iter = driftVolumeMap.find(LArPandoraGeometry::GetTpcID(cstat, tpc));

    if (driftVolumeMap.end() == iter)
        throw cet::exception("LArPandora") << " LArPandoraGeometry::GetVolumeID --- found a TPC that doesn't belong to a drift volume";

    return iter->second.GetVolumeID();
}

//------------------------------------------------------------------------------------------------------------------------------------------

geo::View_t LArPandoraGeometry::GetGlobalView(const unsigned int cstat, const unsigned int tpc, const geo::View_t hit_View)
{
    const bool switchUV(LArPandoraGeometry::ShouldSwitchUV(cstat, tpc));

    if (hit_View == geo::kW)
    {
        return geo::kW;
    }
    else if(hit_View == geo::kU)
    {
        return (switchUV ? geo::kV : geo::kU);
    }
    else if(hit_View == geo::kV)
    {
        return (switchUV ? geo::kU : geo::kV);
    }
    else
    {
        throw cet::exception("LArPandora") << " LArPandoraGeometry::GetGlobalView --- found an unknown plane view (not U, V or W) ";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPandoraGeometry::GetTpcID(const unsigned int cstat, const unsigned int tpc)
{
    // We assume there will never be more than 10000 TPCs in a cryostat!
    if (tpc >= 10000)
        throw cet::exception("LArPandora") << " LArPandoraGeometry::GetTpcID --- found a TPC with an ID greater than 10000 ";

    return 10000 * cstat + tpc;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraGeometry::ShouldSwitchUV(const unsigned int cstat, const unsigned int tpc)
{
    // We determine whether U and V views should be switched by checking the drift direction
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::TPCGeo &theTpc(theGeometry->TPC(tpc, cstat));
    const bool isPositiveDrift(theTpc.DriftDirection() == geo::kPosX);
    return LArPandoraGeometry::ShouldSwitchUV(isPositiveDrift);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraGeometry::ShouldSwitchUV(const bool isPositiveDrift)
{
    // We assume that all multiple drift volume detectors have the APA - CPA - APA - CPA design
    return isPositiveDrift;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::LoadGeometry(LArDriftVolumeList &driftVolumeList)
{
    // This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction,
    // common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
    //
    // ATTN: we assume that all two-wire detectors use the U and V views, and that the wires in the W view are always vertical
    //       An exception will be thrown if W wires are not vertical, indicating that LArPandoraGeometry should not be used.

    if (!driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- detector geometry has already been loaded ";

    typedef std::set<unsigned int> UIntSet;

    // Load Geometry Service
    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int wirePlanes(theGeometry->MaxPlanes());

    const double wirePitchU(theGeometry->WirePitch(geo::kU));
    const double wirePitchV(theGeometry->WirePitch(geo::kV));
    const double wirePitchW((wirePlanes > 2) ? theGeometry->WirePitch(geo::kW) : 0.5 * (wirePitchU + wirePitchV));

    const double maxDeltaTheta(0.01); // leave this hard-coded for now

    // Loop over cryostats
    for (unsigned int icstat = 0; icstat < theGeometry->Ncryostats(); ++icstat)
    {
        UIntSet cstatList;

        // Loop over TPCs in in this cryostat
        for (unsigned int itpc1 = 0; itpc1 < theGeometry->NTPC(icstat); ++itpc1)
        {
            if (cstatList.end() != cstatList.find(itpc1))
                continue;

            // Use this TPC to seed a drift volume
            const geo::TPCGeo &theTpc1(theGeometry->TPC(itpc1, icstat));
            cstatList.insert(itpc1);

            const double wireAngleU(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat));
            const double wireAngleV((0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat)) * -1.f);
            const double wireAngleW((wirePlanes > 2) ? (0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat)) : 0.0);

            if (std::fabs(wireAngleW) > maxDeltaTheta)
                throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- the W-wires are not vertical in this detector ";

            double localCoord1[3] = {0.,0.,0.};
            double worldCoord1[3] = {0.,0.,0.};
            theTpc1.LocalToWorld(localCoord1, worldCoord1);

            const double min1(worldCoord1[0] - 0.5 * theTpc1.ActiveHalfWidth());
            const double max1(worldCoord1[0] + 0.5 * theTpc1.ActiveHalfWidth());

            double driftMinX(worldCoord1[0] - theTpc1.ActiveHalfWidth());
            double driftMaxX(worldCoord1[0] + theTpc1.ActiveHalfWidth());
            double driftMinY(worldCoord1[1] - theTpc1.ActiveHalfHeight());
            double driftMaxY(worldCoord1[1] + theTpc1.ActiveHalfHeight());
            double driftMinZ(worldCoord1[2] - 0.5 * theTpc1.ActiveLength());
            double driftMaxZ(worldCoord1[2] + 0.5 * theTpc1.ActiveLength());

            const bool isPositiveDrift(theTpc1.DriftDirection() == geo::kPosX);

            UIntSet tpcList;
            tpcList.insert(itpc1);

            // Now identify the other TPCs associated with this drift volume
            for (unsigned int itpc2 = itpc1+1; itpc2 < theGeometry->NTPC(icstat); ++itpc2)
            {
                if (cstatList.end() != cstatList.find(itpc2))
                    continue;

                const geo::TPCGeo &theTpc2(theGeometry->TPC(itpc2, icstat));

                if (theTpc1.DriftDirection() != theTpc2.DriftDirection())
                    continue;

                const double dThetaU(theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kU, itpc2, icstat));
                const double dThetaV(theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kV, itpc2, icstat));
                const double dThetaW((wirePlanes > 2) ? (theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kW, itpc2, icstat)) : 0.0);

                if (dThetaU > maxDeltaTheta || dThetaV > maxDeltaTheta || dThetaW > maxDeltaTheta)
                    continue;

                double localCoord2[3] = {0.,0.,0.};
                double worldCoord2[3] = {0.,0.,0.};
                theTpc2.LocalToWorld(localCoord2, worldCoord2);
                const double min2(worldCoord2[0] - 0.5 * theTpc2.ActiveHalfWidth());
                const double max2(worldCoord2[0] + 0.5 * theTpc2.ActiveHalfWidth());

                if ((min2 > max1) || (min1 > max2))
                    continue;

                cstatList.insert(itpc2);
                tpcList.insert(itpc2);

                driftMinX = std::min(driftMinX, worldCoord2[0] - theTpc2.ActiveHalfWidth());
                driftMaxX = std::max(driftMaxX, worldCoord2[0] + theTpc2.ActiveHalfWidth());
                driftMinY = std::min(driftMinY, worldCoord2[1] - theTpc2.ActiveHalfHeight());
                driftMaxY = std::max(driftMaxY, worldCoord2[1] + theTpc2.ActiveHalfHeight());
                driftMinZ = std::min(driftMinZ, worldCoord2[2] - 0.5 * theTpc2.ActiveLength());
                driftMaxZ = std::max(driftMaxZ, worldCoord2[2] + 0.5 * theTpc2.ActiveLength());
            }

            // Collate the tpc volumes in this drift volume
            LArTpcVolumeList tpcVolumeList;

            for(const unsigned int itpc : tpcList)
            {
                tpcVolumeList.push_back(LArTpcVolume(icstat, itpc));
            }

            // Create new daughter drift volume (volume ID = 0 to N-1)
            driftVolumeList.push_back(LArDriftVolume(driftVolumeList.size(), isPositiveDrift,
                wirePitchU, wirePitchV, wirePitchW, wireAngleU, wireAngleV,
                0.5 * (driftMaxX + driftMinX), 0.5 * (driftMaxY + driftMinY), 0.5 * (driftMaxZ + driftMinZ),
                (driftMaxX - driftMinX), (driftMaxY - driftMinY), (driftMaxZ - driftMinZ),
                (wirePitchU + wirePitchV + wirePitchW + 0.1), tpcVolumeList));
        }
    }

    if (driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- failed to find any drift volumes in this detector geometry ";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::LoadGlobalParentGeometry(const LArDriftVolumeList &driftVolumeList, LArDriftVolumeList &parentVolumeList)
{
    // This method will create one or more master drift volumes (these are groups of drift volumes that share a common drift orientation
    // along the X-axis, have parallel or near-parallel wire angles, and similar wire pitches

    if (!parentVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalParentGeometry --- parent geometry has already been loaded ";

    if (driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalParentGeometry --- detector geometry has not yet been loaded ";

    // Take most properties from first volume, then loop over all volumes
    const LArDriftVolume &frontVolume(driftVolumeList.front());

    const bool   switchViews(LArPandoraGeometry::ShouldSwitchUV(frontVolume.IsPositiveDrift()));
    const bool   isPositiveDrift(switchViews ? !frontVolume.IsPositiveDrift() : frontVolume.IsPositiveDrift());
    const double parentWirePitchU(std::max(frontVolume.GetWirePitchU(), frontVolume.GetWirePitchV()));
    const double parentWirePitchV(parentWirePitchU);
    const double parentWirePitchW(frontVolume.GetWirePitchW());
    const double parentWireAngleU(switchViews ? - frontVolume.GetWireAngleV() : frontVolume.GetWireAngleU());
    const double parentWireAngleV(switchViews ? - frontVolume.GetWireAngleU() : frontVolume.GetWireAngleV());

    double parentMinX(frontVolume.GetCenterX() - 0.5 * frontVolume.GetWidthX());
    double parentMaxX(frontVolume.GetCenterX() + 0.5 * frontVolume.GetWidthX());
    double parentMinY(frontVolume.GetCenterY() - 0.5 * frontVolume.GetWidthY());
    double parentMaxY(frontVolume.GetCenterY() + 0.5 * frontVolume.GetWidthY());
    double parentMinZ(frontVolume.GetCenterZ() - 0.5 * frontVolume.GetWidthZ());
    double parentMaxZ(frontVolume.GetCenterZ() + 0.5 * frontVolume.GetWidthZ());

    LArTpcVolumeList tpcVolumeList;

    for (const LArDriftVolume &driftVolume : driftVolumeList)
    {
        parentMinX = std::min(parentMinX, driftVolume.GetCenterX() - 0.5 * driftVolume.GetWidthX());
        parentMaxX = std::max(parentMaxX, driftVolume.GetCenterX() + 0.5 * driftVolume.GetWidthX());
        parentMinY = std::min(parentMinY, driftVolume.GetCenterY() - 0.5 * driftVolume.GetWidthY());
        parentMaxY = std::max(parentMaxY, driftVolume.GetCenterY() + 0.5 * driftVolume.GetWidthY());
        parentMinZ = std::min(parentMinZ, driftVolume.GetCenterZ() - 0.5 * driftVolume.GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, driftVolume.GetCenterZ() + 0.5 * driftVolume.GetWidthZ());

        tpcVolumeList.insert(tpcVolumeList.end(), driftVolume.GetTpcVolumeList().begin(), driftVolume.GetTpcVolumeList().end());
    }

    // Create parent drift volume (volume ID = 0)
    parentVolumeList.push_back(LArDriftVolume(0, isPositiveDrift,
        parentWirePitchU, parentWirePitchV, parentWirePitchW, parentWireAngleU, parentWireAngleV,
        0.5 * (parentMaxX + parentMinX), 0.5 * (parentMaxY + parentMinY), 0.5 * (parentMaxZ + parentMinZ),
        (parentMaxX - parentMinX), (parentMaxY - parentMinY), (parentMaxZ - parentMinZ),
        frontVolume.GetSigmaUVZ(), tpcVolumeList));

    if (parentVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalParentGeometry --- failed to create parent geometry list ";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::LoadGlobalDaughterGeometry(const LArDriftVolumeList &driftVolumeList, LArDriftVolumeList &daughterVolumeList)
{
    // This method will create one or more daughter volumes (these share a common drift orientation along the X-axis,
    // have parallel or near-parallel wire angles, and similar wire pitches)
    //
    // ATTN: we assume that the U and V planes have equal and opposite wire orientations

    if (!daughterVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- daughter geometry has already been loaded ";

    if (driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- detector geometry has not yet been loaded ";

    // Create daughter drift volumes
    for (const LArDriftVolume &driftVolume : driftVolumeList)
    {
        const bool   switchViews(LArPandoraGeometry::ShouldSwitchUV(driftVolume.IsPositiveDrift()));
        const double daughterWirePitchU(switchViews ? driftVolume.GetWirePitchV() : driftVolume.GetWirePitchU());
        const double daughterWirePitchV(switchViews ? driftVolume.GetWirePitchU() : driftVolume.GetWirePitchV());
        const double daughterWirePitchW(driftVolume.GetWirePitchW());
        const double daughterWireAngleU(switchViews ? - driftVolume.GetWireAngleV() : driftVolume.GetWireAngleU());
        const double daughterWireAngleV(switchViews ? - driftVolume.GetWireAngleU() : driftVolume.GetWireAngleV());

        daughterVolumeList.push_back(LArDriftVolume(driftVolume.GetVolumeID(), driftVolume.IsPositiveDrift(),
            daughterWirePitchU, daughterWirePitchV, daughterWirePitchW, daughterWireAngleU, daughterWireAngleV,
            driftVolume.GetCenterX(), driftVolume.GetCenterY() , driftVolume.GetCenterZ(),
            driftVolume.GetWidthX(), driftVolume.GetWidthY(), driftVolume.GetWidthZ(),
            driftVolume.GetSigmaUVZ(), driftVolume.GetTpcVolumeList()));
    }

    if (daughterVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- failed to create daughter geometry list ";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::WritePrecisionElement(pandora::TiXmlElement *const pParentElement, const std::string &elementName, const double value)
{
    const std::streamsize ss(std::cout.precision());
    std::ostringstream oss;

    if ((oss << std::setprecision(12) << value << std::setprecision(ss)).fail())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::WritePrecisionElement --- Could not write LArDriftVolumes to xml file, TypeToString failed.";

    pandora::TiXmlElement *const pDaughterElement = new pandora::TiXmlElement(elementName);
    pDaughterElement->LinkEndChild(new pandora::TiXmlText(oss.str()));
    pParentElement->LinkEndChild(pDaughterElement);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDriftVolume::LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
    const double wirePitchU, const double wirePitchV, const double wirePitchW, const double wireAngleU, const double wireAngleV,
    const double centerX, const double centerY, const double centerZ, const double widthX, const double widthY, const double widthZ,
    const double sigmaUVZ, const LArTpcVolumeList tpcVolumeList) :
    m_volumeID(volumeID),
    m_isPositiveDrift(isPositiveDrift),
    m_wirePitchU(wirePitchU),
    m_wirePitchV(wirePitchV),
    m_wirePitchW(wirePitchW),
    m_wireAngleU(wireAngleU),
    m_wireAngleV(wireAngleV),
    m_centerX(centerX),
    m_centerY(centerY),
    m_centerZ(centerZ),
    m_widthX(widthX),
    m_widthY(widthY),
    m_widthZ(widthZ),
    m_sigmaUVZ(sigmaUVZ),
    m_tpcVolumeList(tpcVolumeList)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArTpcVolumeList &LArDriftVolume::GetTpcVolumeList() const
{
    return m_tpcVolumeList;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraGeometry::Settings::Settings() :
    m_globalCoordinates(false),
    m_globalDriftVolume(false),
    m_printGeometry(false)
{
}

} // namespace lar_pandora
