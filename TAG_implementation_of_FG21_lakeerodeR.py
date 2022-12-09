#!/usr/bin/env python
#TAG implementation of CIF lakeerodeR.py code. New code added by TAG will be explicitly commented as added by TAG
import os, optparse, subprocess, sys, tempfile, ogr, osr, gdal

import copy
import anuga
import numpy as np
import math
from suspendedtransport_operator import suspendedtransport_operator
from bedloadtransport_operator import bedloadtransport_operator
from friction_operator import friction_operator
from AngleofRepose_operator import AoR_operator
from analysisfunctions import finalstats
from model_params import gravity, grainD, constantn


def man(option, opt, value, parser):
    print >>sys.stderr, parser.usage
    print >>sys.stderr, '''\
This program does lake erosion in ANUGA.  
'''
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
    try:
        try:
            usage = "usage: TAG_implementation_of_FG21_lakeerodeR.py LakeRadius-in-meters LakeDepth-in-meters TopoExponent RimHeight-in-m RegionalSlope-in-unitless [-i initialbreachdepth]\n" #TAG call, which requires lake radius, lake depth, lake topography power law exponent, rim height, regional slope (m/m), and optionally the initial depth of the notch in the rim (i.e., the initial breach depth)
            parser = optparse.OptionParser(usage=usage)
            parser.add_option("-i", dest="initbreachdepth", help="Initial breach depth for flood")
            (options, inargs) = parser.parse_args()
            if not inargs: parser.error("need parameters")
            lakeradius=float(inargs[0])
            h=float(inargs[1]) #TAG - get lake depth
            lakeexponent=float(inargs[2]) #TAG - get lake topopgraphy power law exponent
            rimheight=float(inargs[3]) #TAG - get rim height
            extslope=-1.0*abs(float(inargs[4])) #TAG - get exterior slope

            rimwidth_sloping=0.04*lakeradius #TAG - set downstream 'width' of sloping part of the rim, not dependant on rim height, but is dependant on lake radius.
            rimwidth_flat=0.5*rimwidth_sloping #TAG - set downstream 'width' of flat part of rim, not dependant on rim height, but is dependant on lake radius.
            rimwidth=rimwidth_flat+rimwidth_sloping #TAG - set downstream 'width' of full rim, not dependant on rim height, but is dependant on lake radius.
            rimslope=-1.0*rimheight/rimwidth_sloping #TAG - Calculate the slope of the crater rim.

            initbreachdepth=5.0          #default
            lakebelowconfining=5.0
            if( options.initbreachdepth): initbreachdepth=float(options.initbreachdepth)
            initbreachdepth=initbreachdepth+lakebelowconfining                            #set this to be a real initial head, not dependent on lakebelowconfining 
            extslopeint=0.01 #TAG - set cross-valley slope to be the same as the regional slope in the base case = 0.01 [-]
 
            #anuga critical parameter
            anuga.g=gravity


            #TAG - calculate length of downstream receiving basin (flat area). Want length such that this lake would only raise 1 m in water level if the full crater lake drained into it.
            lakemaxcontour=-1.0*lakebelowconfining
            nominal_lake_vol=2.0*math.pi*(((lakemaxcontour+h)/2.0)*((h+lakemaxcontour)*(lakeradius**lakeexponent)/h)**(2.0/lakeexponent) - (h/((2.0+lakeexponent)*lakeradius**lakeexponent))*((h+lakemaxcontour)*(lakeradius**lakeexponent)/h)**((2.0+lakeexponent)/lakeexponent))
            dslakelength=nominal_lake_vol/(2.0*lakeradius)
            print "NOMINAL LAKE VOLUME:"
            print nominal_lake_vol

            #TAG - set the sloping runout length - length is set so that the cross-sectional area is the same for all runs, which ultimately is set so that at the minimum regional slope value the receiving basin is 0.5*h below the lake base (i.e., slope drop = 1.5*h - rimheight). Also want to set it so the high-resolution mesh is only for 3*lake_radius.
            sloperunoutlength=3.0 #TAG - this value is not absolute, but relative to the lake radius
            Smin=0.001
            downstreamsloperunoutlength=((((((1.5*h-rimheight)**2)/(abs(extslope)*Smin))**(0.5)))/lakeradius)-3.0 #TAG - this value is not absolute, but relative to the lake radius
            if (downstreamsloperunoutlength<0): downstreamsloperunoutlength = 0.0 #TAG - Can't be negative., so just means the full runout region will be high-resolution
            fullsloperunoutlength=sloperunoutlength+downstreamsloperunoutlength #TAG - this value is not absolute, but relative to the lake radius

 
            #domain parameters
            domainpolygon=[[0,0],[0,2*lakeradius],[2*lakeradius,2*lakeradius],[2*lakeradius,1.2*lakeradius],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius),1.2*lakeradius],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius),2*lakeradius],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius+dslakelength),2*lakeradius],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius+dslakelength),0],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius),0],[(2*lakeradius+rimwidth+fullsloperunoutlength*lakeradius),0.8*lakeradius],[2*lakeradius,0.8*lakeradius],[2*lakeradius,0],[0,0]]

            boundary_tags={'exterior': [0]}

            high_resolution_init=(lakeradius/100)**2 #TAG - Initial high-resolution grid cell region, upon which other calculations are based.
            base_resolution=high_resolution_init*80.0
            high_resolution=high_resolution_init/2.0
            mid_resolution=high_resolution_init/0.5 #TAG - Medium-resolution grid cell region.
            initbreachwidth=lakeradius/50.0

            #TAG - Set mesh, but need to know if a medium-resolution grid cell region is required, or if the whole runout region is high-resolution.
            if (downstreamsloperunoutlength>0):
                hiresregion=[[1.5*lakeradius,1.2*lakeradius],[1.5*lakeradius,0.8*lakeradius],[(2.0+sloperunoutlength)*lakeradius+rimwidth,0.8*lakeradius],[(2.0+sloperunoutlength)*lakeradius+rimwidth,1.2*lakeradius],[1.5*lakeradius,1.2*lakeradius]]
                midresregion=[[(2.0+sloperunoutlength)*lakeradius+rimwidth,0.8*lakeradius],[(2.0+sloperunoutlength)*lakeradius+rimwidth,1.2*lakeradius],[(2.5+fullsloperunoutlength)*lakeradius+rimwidth,1.2*lakeradius],[(2.5+fullsloperunoutlength)*lakeradius+rimwidth,0.8*lakeradius],[(2.0+sloperunoutlength)*lakeradius+rimwidth,0.8*lakeradius]]
                interior_regions = [[hiresregion, high_resolution],[midresregion,mid_resolution]]
            if (downstreamsloperunoutlength<=0):
                hiresregion=[[1.5*lakeradius,1.2*lakeradius],[1.5*lakeradius,0.8*lakeradius],[(2.5+sloperunoutlength)*lakeradius+rimwidth,0.8*lakeradius],[(2.5+sloperunoutlength)*lakeradius+rimwidth,1.2*lakeradius],[1.5*lakeradius,1.2*lakeradius]]
                interior_regions = [[hiresregion, high_resolution]]


            meshname='lake.msh'
            m = anuga.create_mesh_from_regions(domainpolygon, boundary_tags, maximum_triangle_area=base_resolution, interior_regions=interior_regions, filename=meshname, use_cache=False)
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration','sedvolx','sedvoly']
            domain = anuga.Domain(meshname, use_cache=False, evolved_quantities = evolved_quantities)
            domain.g = anuga.g  # make sure the domain inherits the package's gravity.

            print 'Number of triangles = ', len(domain)

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return idx
                
            def topography(x, y):
                xc=x-lakeradius 
                yc=y-lakeradius

                #TAG - Setup lake topography according to a power law.
                lakepowerlaw=h*(((((xc**2+yc**2)**0.5)/lakeradius)**lakeexponent)-1)
                z=np.where(lakepowerlaw<0,lakepowerlaw,0.0)

                #TAG - Setup flat portion of the rim.
                z=np.where(np.logical_and(x>lakeradius*2,x<=(lakeradius*2+rimwidth_flat)),abs(y-lakeradius)*extslopeint,z)

                #TAG - Setup steep portion of the rim.
                z=np.where(np.logical_and(x>(lakeradius*2+rimwidth_flat),x<=(lakeradius*2+rimwidth)),(x-lakeradius*2-rimwidth_flat)*rimslope+abs(y-lakeradius)*extslopeint,z)

                #TAG - Setup backround slope.
                z=np.where(np.logical_and(x>(lakeradius*2+rimwidth),x<=(lakeradius*(2+fullsloperunoutlength)+rimwidth)),(x-lakeradius*2-rimwidth)*extslope+abs(y-lakeradius)*extslopeint-abs(np.amin(np.where(np.logical_and(x>lakeradius*2,x<=(lakeradius*2+rimwidth)),z,0.0))),z)

                #TAG - Setup receiving basin.
                z=np.where(x>(lakeradius*(2+fullsloperunoutlength)+rimwidth),np.amin(np.where(np.logical_and(x>(lakeradius*2+rimwidth),x<=(lakeradius*(2+fullsloperunoutlength)+rimwidth)),z,0.0))+abs(y-lakeradius)*extslopeint,z)

                # add gaussian noise, with zero mean and 0.1 m standard deviation.
                #z=z+np.random.normal(0,0.1,z.shape)
                
                # set up breach as a failure at the midpoint in a region
                #mask1 = (x>=(lakeradius*2+(initbreachdepth/extslope)))
                mask1 = (x>=lakeradius)
                #mask2 = (x<=(lakeradius*2-(initbreachdepth/extslope)))
                mask3 = (y>=((lakeradius-(initbreachwidth/2.0))))  
                mask4 = (y<=((lakeradius+(initbreachwidth/2.0))))
                mask5 = (z>-initbreachdepth)
                #mask=mask1 & mask2 & mask3 & mask4 & mask5
                mask=mask1 & mask3 & mask4 & mask5
                z[mask]=-initbreachdepth


                return z

            def initialstage(x,y):
                clow=domain.get_quantity('elevation').get_values(interpolation_points=[[2*lakeradius,lakeradius]])
                istage=clow+initbreachdepth-lakebelowconfining
                
                #outsidelake. Force depth=0 (istage=topo) across rim.
                #istage=np.where(x>lakeradius*2,topography(x,y),istage)
                istage=np.where(x>lakeradius*2,-10000,istage)

                return istage

            name="Marslake"+str(lakeradius)+"_"+str(grainD*1000)+"mm"
            domain.set_name(name)
            # Choose between DE0 (less accurate), DE1 (more accurate, slower), DE2 (even more accurate, slower still)
            domain.set_flow_algorithm('DE0')
            domain.set_CFL(cfl=1)
            domain.set_minimum_allowed_height(0.01)  #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for computation)
            domain.set_minimum_storable_height(0.01) #raise the default minimum depth from 1 mm to 1 cm to speed the simulation (for storage/visualization)
            domain.set_maximum_allowed_speed(1)      #maximum particle speed that is allowed in water shallower than minimum_allowed_height (back of the envelope suggests 1m/s should be below tau_crit)

 

            np.seterr(invalid='ignore')
            domain.set_quantity('elevation', topography, location='centroids') # elevation is a function
            voltotal=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
            domain.set_quantity('concentration', 0.00)        # Start with no esd
            domain.set_quantity('friction', 0.0545)      # Constant Manning friction.  Only used to initialize domain (calculated every time step).
            domain.set_quantity('stage', initialstage, location='centroids')
            domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 2})

            np.seterr(invalid='warn')

            # Setup boundary conditions

            Br=anuga.Reflective_boundary(domain)
            domain.set_boundary({'exterior': Br})

            operatezero = suspendedtransport_operator(domain)
            operateone = bedloadtransport_operator(domain)
            operatetwo = friction_operator(domain)
            operatethree = AoR_operator(domain)

            initial=domain.get_quantity('elevation').get_values(interpolation_points=[[2.01*lakeradius,lakeradius]], location='centroids')
            initialdomain=copy.deepcopy(domain)
            ystep=300
            ftime=86400*300
            count=0
            totdv=0.00001 #TAG - This is the starting value for the total volume drained, ultimately used as a threshold for cutoff. Need to set a non-zero value, or it throws an error.
            estinstdv_thresh=0.00001 #TAG - This is the threshold value for the volume drained in one time step, ultimately used as a cutoff.
            estinstdv_thresh_count=0 #TAG - This is the counter for number of continuous time steps the drained volume exceeds the threshold value.
            estinstdv_thresh_count_stop=12 #TAG - After how many continuous timesteps below estinstdv_thresh do you shutdown; 12 = 1 h
            erosion_cutoff=0.001 #TAG - This is the threshold value for the erosion cutoff = 1 mm/timestep.
            erosion_cutoff_start=0.05 #TAG - This defines when to start monitoring for the erosion cutoff (since it takes some time for the flood to get going, and don't want to cutoff then) = 5 cm/timestep.
            last_step_sill_elev=-1.0*initbreachdepth #TAG - Start the breach floor elevation at init breach depth.
            timestep_erosion=0.0 #TAG - This is how much erosion at the breach occured during this timestep. Start at 0.0.
            erosion_cutoff_counter=0 #TAG - This is the counter for number of continuous time steps the erosion rate exceeds the threshold value.
            erosion_cutoff_counter_stop=288 #TAG - After how many continuous timesteps below erosion_cutoff do you shutdown; 288 = 24 h
            ec_t_threshold=0 #TAG - Flag to indicate that you don't start evaluating the erosion cutoff until it goes above erosion_cutoff_start the first time
            #TAG - Setup a polyline halfway into the flat part of the crater rim, which is = 2.0*lake_radius+0.5*rimwidth_flat, where rimwidth_flat=0.5*rimwidth_sloping and rimwidth_sloping=0.04*lakeradius, so 2.01*lakeradius. Also use 0.8/1.2*lake_radius as the y points, since this is what is used in the analysis function (this is the narrow bit of the domain).
            acrossthebreach=[[2.01*lakeradius,0.8*lakeradius],[2.01*lakeradius,1.2*lakeradius]]
            #TAG - Set the sill elevation profile spacing (both x=along profile and y=between profiles) as the high-res mesh size triangle equilateral-equivalent side length, divided by 100 - massive oversampling, but that is OK :). Same as across-breach profile spacing in analysis function.
            sill_profile_crossspacing=(((((lakeradius/100)**2)/2.0)*(4/(3**(0.5))))**(0.5))/10.0
            #TAG - Set x locations of sill profile --> goes from lake center to halfway through the flat part of the crater rim.
            sill_profile_x_locs=np.arange(lakeradius,2.01*lakeradius,sill_profile_crossspacing)
            #TAG - Set y locations of sill profile
            sill_profile_y_locs=np.arange((lakeradius-sill_profile_crossspacing),(lakeradius+2*sill_profile_crossspacing+1),sill_profile_crossspacing)
            #TAG - Set xy profile pairs
            sill_profile_xy_pairs=np.column_stack((sill_profile_x_locs,(np.zeros(len(sill_profile_x_locs))+(lakeradius-2*sill_profile_crossspacing))))
            for temp_y in sill_profile_y_locs:
                sill_profile_xy_pairs=np.vstack((sill_profile_xy_pairs,np.column_stack((sill_profile_x_locs,(np.zeros(len(sill_profile_x_locs))+temp_y)))))
            t_threshold=20000 #TAG - Let run go at LEAST this long = ~5.5 hrs.

            for t in domain.evolve(yieldstep=ystep, finaltime=ftime):
                print domain.timestepping_statistics()
                print 'xmom:'+str(domain.get_quantity('xmomentum').get_values(interpolation_points=[[2*lakeradius,lakeradius]], location='centroids'))   
                volcurr=np.sum(domain.quantities['elevation'].centroid_values*domain.areas)
                #TAG get sill elevation
                sillprofs=domain.get_quantity('elevation').get_values(interpolation_points=sill_profile_xy_pairs)
                sillelev=np.amax(sillprofs)
                timestep_erosion= last_step_sill_elev - sillelev #TAG - Should always be going down, so this should always be positive...
                last_step_sill_elev=sillelev #TAG - Swap over last sill elev
                print 'breach sill elevation: '+str(sillelev)
                print 'breach sill dz/dt: '+str(timestep_erosion)
                volsed=np.sum(domain.quantities['concentration'].centroid_values*domain.quantities['height'].centroid_values*domain.areas)            
                conservation=(volcurr+volsed-voltotal)/voltotal
                print 'conservation: '+'{:.8%}'.format(conservation)
                breachcrosssectionflow=abs(domain.get_flow_through_cross_section(acrossthebreach))
                estinstdv=breachcrosssectionflow*ystep
                totdv=totdv+estinstdv
                estinstdv_frac=estinstdv/totdv
                print 'Q: '+str(breachcrosssectionflow)
                print 'Inst_DV_Fraction: '+str(estinstdv_frac)

                if (timestep_erosion>erosion_cutoff_start)&(ec_t_threshold==0)&(t>t_threshold):
                    #TAG - Go above erosion cutoff start threshold first time (after t_threshold), i.e., runaway erosion has likely started, so can start monitoring for shutdown conditions.
                    ec_t_threshold=1
                    print 'Passed sill erosion threshold. Staring to look for shutoff.'

                if (estinstdv_frac<estinstdv_thresh)&(ec_t_threshold==1):
                    estinstdv_thresh_count=estinstdv_thresh_count+1
                    print 'Dropped below the instantaneous theshold drained volume, time '+str(estinstdv_thresh_count)+' of '+str(estinstdv_thresh_count_stop)
                elif (estinstdv_thresh_count>0)&(ec_t_threshold==1)&(estinstdv_frac>estinstdv_thresh):
                    #TAG - NOT below threshold, but count has been set before, so reset to 0.
                    estinstdv_thresh_count=0
                
                if (timestep_erosion<erosion_cutoff)&(ec_t_threshold==1):
                    erosion_cutoff_counter=erosion_cutoff_counter+1
                    print 'Dropped below the threshold sill dz/dt, time '+str(erosion_cutoff_counter)+' of '+str(erosion_cutoff_counter_stop)
                elif (timestep_erosion>erosion_cutoff)&(ec_t_threshold==1)&(erosion_cutoff_counter>0):
                    #TAG - NOT below erosion dz/dt threshold, but count has been set before, so reset to 0.
                    erosion_cutoff_counter=0

                if (erosion_cutoff_counter==erosion_cutoff_counter_stop)&(t>t_threshold):
                    print "Dropped below the sill dz/dt threshold for the final time. Ending run."
                    break

                if (estinstdv_thresh_count==estinstdv_thresh_count_stop)&(t>t_threshold):
                    print "Dropped below the instantaneous theshold drained volume for the final time. Ending run."
                    break
                count=count+1

            time, Q = anuga.get_flow_through_cross_section(name+'.sww',acrossthebreach)
            print Q

            initname="initial_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".asc"
            finname="final_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".asc"
            fluxname="flux_"+str(lakeradius)+"_"+str(grainD*1000)+"mm"+".txt"
            np.savetxt(fluxname,Q)
            finalstats(domain,initialdomain,lakeradius)
            np.save('XconcC.npy', domain.quantities['concentration'].centroid_values)
            np.save('XelevC.npy', domain.quantities['elevation'].centroid_values)
            np.save('XxmC.npy', domain.quantities['xmomentum'].centroid_values)
            np.save('XymC.npy', domain.quantities['ymomentum'].centroid_values)
            np.save('XstageC.npy', domain.quantities['stage'].centroid_values)
            np.save('XconcV.npy', domain.quantities['concentration'].vertex_values)
            np.save('XelevV.npy', domain.quantities['elevation'].vertex_values)
            np.save('XxmV.npy', domain.quantities['xmomentum'].vertex_values)
            np.save('XymV.npy', domain.quantities['ymomentum'].vertex_values)
            np.save('XstageV.npy', domain.quantities['stage'].vertex_values)

        except optparse.OptionError, msg:
            raise Usage(msg)

    except Usage, err:
        print >>sys.stderr, err.msg
        # print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())
