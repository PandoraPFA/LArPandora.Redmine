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

namespace lar_pandora
{

void LArPandoraGeometry::LoadGeometry(LArDriftVolumeList &driftVolumeList)
{
    // This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction, 
    // common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
    // Note: we assume that all two-wire detectors use the U and V views, and that the wires in the W view are always vertical

    if (!driftVolumeList.empty())
        throw cet::exception("LArPandora") << " Throwing exception - detector geoemetry has already been loaded ";

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
                throw cet::exception("LArPandora") << " Throwing exception - the W-wires are not vertical in this detector ";

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
            for(UIntSet::const_iterator iter = tpcList.begin(), iterEnd = tpcList.end(); iter != iterEnd; ++iter)
	    {
	        const unsigned int itpc = *iter;
	        tpcVolumeList.push_back(LArTpcVolume(icstat, itpc));
	    }

	    // Create the new drift volume
            driftVolumeList.push_back(LArDriftVolume(driftVolumeList.size(), isPositiveDrift,
                wirePitchU, wirePitchV, wirePitchW, wireAngleU, wireAngleV,
	        0.5 * (driftMaxX + driftMinX), 0.5 * (driftMaxY + driftMinY), 0.5 * (driftMaxZ + driftMinZ),
		(driftMaxX - driftMinX), (driftMaxY - driftMinY), (driftMaxZ - driftMinZ),
		(wirePitchU + wirePitchV + wirePitchW + 0.1), tpcVolumeList));
	}
    }

    if (driftVolumeList.empty())
        throw cet::exception("LArPandora") << " Throwing exception - failed to find any drift volumes in this detector geometry ";
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraGeometry::PrintGeometry(const LArDriftVolumeList &driftVolumeList)
{
    // Useful for debugging (or reverse-engineering) a detector geometry
    std::cout << " *** LArPandoraModule::PrintGeometry() *** " << std::endl;

    for (LArDriftVolumeList::const_iterator iter1 = driftVolumeList.begin(), iterEnd1 = driftVolumeList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArDriftVolume &driftVolume = *iter1;

	std::cout << " *** DriftVolume *** " << std::endl;
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

        const LArTpcVolumeList &tpcVolumeList = driftVolume.GetTpcVolumeList();

        for (LArTpcVolumeList::const_iterator iter2 = tpcVolumeList.begin(), iterEnd2 = tpcVolumeList.end(); iter2 != iterEnd2; ++iter2)
	{
	    const LArTpcVolume &tpcVolume = *iter2; 
	    std::cout << "   [" << tpcVolume.GetCryostat() << "][" << tpcVolume.GetTpc() << "]" << std::endl;
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDriftVolume::LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
  const double wirePitchU, const double wirePitchV, const double wirePitchW, const double wireAngleU, const double wireAngleV,
  const double centerX, const double centerY, const double centerZ, const double widthX, const double widthY, const double widthZ,
  const double sigmaUVZ, const LArTpcVolumeList tpcVolumeList) :
    m_volumeID(volumeID), m_isPositiveDrift(isPositiveDrift), m_wirePitchU(wirePitchU), m_wirePitchV(wirePitchV), m_wirePitchW(wirePitchW),
    m_wireAngleU(wireAngleU), m_wireAngleV(wireAngleV), m_centerX(centerX), m_centerY(centerY), m_centerZ(centerZ),
    m_widthX(widthX), m_widthY(widthY), m_widthZ(widthZ), m_sigmaUVZ(sigmaUVZ), m_tpcVolumeList(tpcVolumeList)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArTpcVolumeList &LArDriftVolume::GetTpcVolumeList() const
{
    return m_tpcVolumeList;
}

} // namespace lar_pandora
