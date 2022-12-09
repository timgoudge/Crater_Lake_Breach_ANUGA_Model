""" -----------------------------------------------
Analysis functions
Fassett implementation
======= """
import numpy as np
import anuga



def finalstats(domain, initialdomain, lakeradius, outfile='outstats.csv'):
    crossspacing=(((((lakeradius/100)**2)/2.0)*(4/(3**(0.5))))**(0.5))/10.0   #Set the profile spacing as the high-resolution mesh size triangle equilateral-equivalent side length, divided by 100 - this is a massive oversampling, but that is OK since it will give the same results.
    breach_profile_x_loc=2.01*lakeradius #Set the breach profile location to be halfway into the flat part of the rim, which is where Qw is calculated for the threshold cutoff. So, breach_profile_x_loc = 2.0*lake_radius+0.5*rimwidth_flat, and rimwidth_flat=0.5*rimwidth_sloping, where rimwidth_sloping=0.04*lakeradius. Thus breach_profile_x_loc=2.01*lakeradius.
    acrossthebreach=np.arange(0.8*lakeradius,1.2*lakeradius,crossspacing)
    breacherosion=0
    areasum=0
    for y in acrossthebreach:
        ybreacherosion=domain.get_quantity('elevation').get_values(interpolation_points=[[breach_profile_x_loc,y]], location='centroids')-initialdomain.get_quantity('elevation').get_values(interpolation_points=[[breach_profile_x_loc,y]], location='centroids')
        if ybreacherosion<0: areasum=areasum+float(ybreacherosion)*crossspacing    #rectangular integration  :)
        if ybreacherosion<breacherosion: breacherosion=float(ybreacherosion)
    
    areasum=-areasum  #because negative areas seem weird
    
    elevdiff=domain.quantities['elevation'].centroid_values-initialdomain.quantities['elevation'].centroid_values  #negative values are eroded, positive values are deposited
    finalremovedvolume=-np.sum(elevdiff[elevdiff<0]*domain.areas[elevdiff<0])
    
    transportedmask=domain.get_centroid_coordinates()[:,0]>2*lakeradius
    heights=(domain.quantities['stage'].centroid_values-domain.quantities['elevation'].centroid_values)
    transportedvolume=(np.where(heights>0,heights,0)*domain.areas)[transportedmask].sum()

    tag_stage_diff=initialdomain.quantities['stage'].centroid_values-domain.quantities['stage'].centroid_values  #positive values are stage drop, which should only be in the lake
    tag_drained_vol=np.sum(tag_stage_diff[tag_stage_diff>0]*domain.areas[tag_stage_diff>0])

    
    print 'Breach Erosion Depth (m): '+str(-breacherosion)
    print 'Cross-sectional Area (m^2): '+str(areasum)
    print 'Outlet volume (m^3): '+str(finalremovedvolume)
    print 'Drained volume (m^3):' +str(transportedvolume)
    print 'TAG Drained volume (m^3):' +str(tag_drained_vol)
    with open(outfile, 'a') as the_file:
        the_file.write(domain.get_name()+','+str(-breacherosion)+','+str(areasum)+','+str(finalremovedvolume)+','+str(transportedvolume)+'\n')
