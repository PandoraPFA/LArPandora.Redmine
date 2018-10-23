/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.cxx
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include <iomanip>
#include <set>
#include <sstream>

namespace lar_pandora
{

void LArPandoraGeometry::LoadDetectorGaps(LArDetectorGapList &listOfGaps)
{
    // Detector gaps can only be loaded once - throw an exception if the output lists are already filled
    if (!listOfGaps.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadDetectorGaps --- the list of gaps already exists ";

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

            const float maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?
            const float deltaZ(std::fabs(driftVolume1.GetCenterZ() - driftVolume2.GetCenterZ()));
            const float deltaY(std::fabs(driftVolume1.GetCenterY() - driftVolume2.GetCenterY()));
            const float deltaX(std::fabs(driftVolume1.GetCenterX() - driftVolume2.GetCenterX()));
            const float widthX(0.5f * (driftVolume1.GetWidthX() + driftVolume2.GetWidthX()));
            const float gapX(deltaX - widthX);

            if (gapX < 0.f || gapX > maxDisplacement || deltaY > maxDisplacement || deltaZ > maxDisplacement)
                continue;

            const float X1((driftVolume1.GetCenterX() < driftVolume2.GetCenterX()) ? (driftVolume1.GetCenterX() + 0.5f * driftVolume1.GetWidthX()) :
                (driftVolume2.GetCenterX() + 0.5f * driftVolume2.GetWidthX()));
            const float X2((driftVolume1.GetCenterX() > driftVolume2.GetCenterX()) ? (driftVolume1.GetCenterX() - 0.5f * driftVolume1.GetWidthX()) :
                (driftVolume2.GetCenterX() - 0.5f * driftVolume2.GetWidthX()));
            const float Y1(std::min((driftVolume1.GetCenterY() - 0.5f * driftVolume1.GetWidthY()),
                (driftVolume2.GetCenterY() - 0.5f * driftVolume2.GetWidthY())));
            const float Y2(std::max((driftVolume1.GetCenterY() + 0.5f * driftVolume1.GetWidthY()),
                (driftVolume2.GetCenterY() + 0.5f * driftVolume2.GetWidthY())));
            const float Z1(std::min((driftVolume1.GetCenterZ() - 0.5f * driftVolume1.GetWidthZ()),
                (driftVolume2.GetCenterZ() - 0.5f * driftVolume2.GetWidthZ())));
            const float Z2(std::max((driftVolume1.GetCenterZ() + 0.5f * driftVolume1.GetWidthZ()),
                (driftVolume2.GetCenterZ() + 0.5f * driftVolume2.GetWidthZ())));

            listOfGaps.push_back(LArDetectorGap(X1, Y1, Z1, X2, Y2, Z2));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::LoadGeometry(LArDriftVolumeList &outputVolumeList, LArDriftVolumeMap &outputVolumeMap)
{
    if (!outputVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- the list of drift volumes already exists ";

    // Use a global coordinate system but keep drift volumes separate
    LArDriftVolumeList inputVolumeList;
    LArPandoraGeometry::LoadGeometry(inputVolumeList);
    LArPandoraGeometry::LoadGlobalDaughterGeometry(inputVolumeList, outputVolumeList);

    // Create mapping between tpc/cstat labels and drift volumes
    for (const LArDriftVolume &driftVolume : outputVolumeList)
    {
        for (const LArDaughterDriftVolume &tpcVolume : driftVolume.GetTpcVolumeList())
        {
            (void) outputVolumeMap.insert(LArDriftVolumeMap::value_type(LArPandoraGeometry::GetTpcID(tpcVolume.GetCryostat(), tpcVolume.GetTpc()), driftVolume));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPandoraGeometry::GetVolumeID(const LArDriftVolumeMap &driftVolumeMap, const unsigned int cstat, const unsigned int tpc)
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

    // ATTN This implicitly assumes that there will be u, v and (maybe) one of either w or y views
    if ((hit_View == geo::kW) || (hit_View == geo::kY))
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

    return ((10000 * cstat) + tpc);
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
    // common range of x coordinates, and common detector parameters such as wire pitch and wire angle).
    if (!driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- detector geometry has already been loaded ";

    typedef std::set<unsigned int> UIntSet;

    // Pandora requires three independent images, and ability to correlate features between images (via wire angles and transformation plugin).
    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int nWirePlanes(theGeometry->MaxPlanes());

    if (nWirePlanes > 3)
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- More than three wire planes present ";

    // We here check the plane information only for the first tpc in the first cryostat.
    if ((0 == theGeometry->Ncryostats()) || (0 == theGeometry->NTPC(0)))
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- unable to access first tpc in first cryostat ";

    std::unordered_set<geo::_plane_proj> planeSet;
    for (unsigned int iPlane = 0; iPlane < nWirePlanes; ++iPlane)
        (void) planeSet.insert(theGeometry->TPC(0, 0).Plane(iPlane).View());

    if ((nWirePlanes != planeSet.size()) || !planeSet.count(geo::kU) || !planeSet.count(geo::kV) || (planeSet.count(geo::kW) && planeSet.count(geo::kY)))
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- expect to find u and v views; if there is one further view, it must be w or y ";

    const bool useYPlane((nWirePlanes > 2) && planeSet.count(geo::kY));

    const float wirePitchU(theGeometry->WirePitch(geo::kU));
    const float wirePitchV(theGeometry->WirePitch(geo::kV));
    const float wirePitchW((nWirePlanes < 3) ? 0.5f * (wirePitchU + wirePitchV) : (useYPlane) ? theGeometry->WirePitch(geo::kY) :
        theGeometry->WirePitch(geo::kW));

    const float maxDeltaTheta(0.01f); // leave this hard-coded for now

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

            const float wireAngleU(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat));
            const float wireAngleV(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat));
            const float wireAngleW((nWirePlanes < 3) ? 0.f : (useYPlane) ? (std::fabs(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kY, itpc1, icstat))) :
                (0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat)));

            double localCoord1[3] = {0., 0., 0.};
            double worldCoord1[3] = {0., 0., 0.};
            theTpc1.LocalToWorld(localCoord1, worldCoord1);

            const double min1(worldCoord1[0] - 0.5 * theTpc1.ActiveHalfWidth());
            const double max1(worldCoord1[0] + 0.5 * theTpc1.ActiveHalfWidth());

            float driftMinX(worldCoord1[0] - theTpc1.ActiveHalfWidth());
            float driftMaxX(worldCoord1[0] + theTpc1.ActiveHalfWidth());
            float driftMinY(worldCoord1[1] - theTpc1.ActiveHalfHeight());
            float driftMaxY(worldCoord1[1] + theTpc1.ActiveHalfHeight());
            float driftMinZ(worldCoord1[2] - 0.5f * theTpc1.ActiveLength());
            float driftMaxZ(worldCoord1[2] + 0.5f * theTpc1.ActiveLength());

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

                const float dThetaU(theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kU, itpc2, icstat));
                const float dThetaV(theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kV, itpc2, icstat));
                const float dThetaW((nWirePlanes < 3) ? 0.f : (useYPlane) ? (theGeometry->WireAngleToVertical(geo::kY, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kY, itpc2, icstat)) :
                    (theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kW, itpc2, icstat)));

                if (dThetaU > maxDeltaTheta || dThetaV > maxDeltaTheta || dThetaW > maxDeltaTheta)
                    continue;

                double localCoord2[3] = {0., 0., 0.};
                double worldCoord2[3] = {0., 0., 0.};
                theTpc2.LocalToWorld(localCoord2, worldCoord2);

                const double min2(worldCoord2[0] - 0.5 * theTpc2.ActiveHalfWidth());
                const double max2(worldCoord2[0] + 0.5 * theTpc2.ActiveHalfWidth());

                if ((min2 > max1) || (min1 > max2))
                    continue;

                cstatList.insert(itpc2);
                tpcList.insert(itpc2);

                driftMinX = std::min(driftMinX, static_cast<float>(worldCoord2[0] - theTpc2.ActiveHalfWidth()));
                driftMaxX = std::max(driftMaxX, static_cast<float>(worldCoord2[0] + theTpc2.ActiveHalfWidth()));
                driftMinY = std::min(driftMinY, static_cast<float>(worldCoord2[1] - theTpc2.ActiveHalfHeight()));
                driftMaxY = std::max(driftMaxY, static_cast<float>(worldCoord2[1] + theTpc2.ActiveHalfHeight()));
                driftMinZ = std::min(driftMinZ, static_cast<float>(worldCoord2[2] - 0.5f * theTpc2.ActiveLength()));
                driftMaxZ = std::max(driftMaxZ, static_cast<float>(worldCoord2[2] + 0.5f * theTpc2.ActiveLength()));
            }

            // Collate the tpc volumes in this drift volume
            LArDaughterDriftVolumeList tpcVolumeList;

            for(const unsigned int itpc : tpcList)
            {
                tpcVolumeList.push_back(LArDaughterDriftVolume(icstat, itpc));
            }

            // Create new daughter drift volume (volume ID = 0 to N-1)
            driftVolumeList.push_back(LArDriftVolume(driftVolumeList.size(), isPositiveDrift,
                wirePitchU, wirePitchV, wirePitchW, wireAngleU, wireAngleV, wireAngleW,
                0.5f * (driftMaxX + driftMinX), 0.5f * (driftMaxY + driftMinY), 0.5f * (driftMaxZ + driftMinZ),
                (driftMaxX - driftMinX), (driftMaxY - driftMinY), (driftMaxZ - driftMinZ),
                (wirePitchU + wirePitchV + wirePitchW + 0.1f), tpcVolumeList));
        }
    }

    if (driftVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGeometry --- failed to find any drift volumes in this detector geometry ";
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
        const float daughterWirePitchU(switchViews ? driftVolume.GetWirePitchV() : driftVolume.GetWirePitchU());
        const float daughterWirePitchV(switchViews ? driftVolume.GetWirePitchU() : driftVolume.GetWirePitchV());
        const float daughterWirePitchW(driftVolume.GetWirePitchW());
        const float daughterWireAngleU(switchViews ? driftVolume.GetWireAngleV() : driftVolume.GetWireAngleU());
        const float daughterWireAngleV(switchViews ? driftVolume.GetWireAngleU() : driftVolume.GetWireAngleV());
        const float daughterWireAngleW(driftVolume.GetWireAngleW());

        daughterVolumeList.push_back(LArDriftVolume(driftVolume.GetVolumeID(), driftVolume.IsPositiveDrift(),
            daughterWirePitchU, daughterWirePitchV, daughterWirePitchW, daughterWireAngleU, daughterWireAngleV, daughterWireAngleW,
            driftVolume.GetCenterX(), driftVolume.GetCenterY() , driftVolume.GetCenterZ(),
            driftVolume.GetWidthX(), driftVolume.GetWidthY(), driftVolume.GetWidthZ(),
            driftVolume.GetSigmaUVZ(), driftVolume.GetTpcVolumeList()));
    }

    if (daughterVolumeList.empty())
        throw cet::exception("LArPandora") << " LArPandoraGeometry::LoadGlobalDaughterGeometry --- failed to create daughter geometry list ";
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDriftVolume::LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
        const float wirePitchU, const float wirePitchV, const float wirePitchW, const float wireAngleU, const float wireAngleV, const float wireAngleW,
        const float centerX, const float centerY, const float centerZ, const float widthX, const float widthY, const float widthZ,
        const float sigmaUVZ, const LArDaughterDriftVolumeList &tpcVolumeList) :
    m_volumeID(volumeID),
    m_isPositiveDrift(isPositiveDrift),
    m_wirePitchU(wirePitchU),
    m_wirePitchV(wirePitchV),
    m_wirePitchW(wirePitchW),
    m_wireAngleU(wireAngleU),
    m_wireAngleV(wireAngleV),
    m_wireAngleW(wireAngleW),
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

const LArDaughterDriftVolumeList &LArDriftVolume::GetTpcVolumeList() const
{
    return m_tpcVolumeList;
}

} // namespace lar_pandora
