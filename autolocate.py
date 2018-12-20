"""
Script to run pseudo-automatic location process that is carried out iteratively
whereby we remove picks with the highest residuals at each iteration until
quality criteria are met:
1) minimum number of P-wave picks
2) minimum number of S-wave picks
and 3) the root-mean-square travel limit

=============================================
Requirements:
    * ObsPy >= 1.1.0
    * ...
=============================================

Te Whare wananga o te upoko o te Ika a Maui
VUW
July 2018
Author: Konstantinos Michailos
"""
from obspy import read_events
import subprocess
import glob
import shutil
import warnings
from obspy.io.nordic.core import write_select
import sys
from obspy.io.nordic.core import read_nordic
from obspy.io.nordic.core import readheader
from obspy.io.nordic.core import nordpick
from obspy.io.nordic.core import readwavename
import numpy as np

def quakeml2sfiles(xml_dir, mseed_dir, sfile_dir):
    """
    Function to match .xml and .mseed files and write out sfiles to a given
    output directory. (.xml and .mseed files should have the same filenames).

    Arguments:
        xml_dir   {str} -- path to the .xml files
        mseed_dir {str} -- path to the .mseed files
        outdir    {str} -- path to write the sfiles
    """

    # List of QuakeML files
    quakeml_list = glob.glob(xml_dir + '2010-04*.xml')
    quakeml_list.sort()
    print('.xml list created')

    # list of mseed files
    mseed_list = glob.glob(mseed_dir + '2010-04*.MSEED')
    mseed_list.sort()
    print('.mseed list created')

    # match times of flist with mlist, remove uneeded waveforms
    # strips file extension and path from filename
    mlist3 = [mlist2.strip('.MSEED') for mlist2 in mseed_list]
    mlist5 = [mlist4.split('/')[-1] for mlist4 in mlist3]
    # the same thing as above
    flist3 = [flist2.strip('.xml') for flist2 in quakeml_list]
    flist5 = [flist4.split('/')[-1] for flist4 in flist3]

    # Create a list of the common elements from the 2 lists (flist5 and mlist5)
    common = []
    for i in range(len(flist5)):
            a = flist5[i]
            for j in range(len(mlist5)):
                b = mlist5[j]
                if a == b:
                    if b not in common:
                        common.append(b)
                    else:
                        continue
                else:
                    continue
    # Remove the elements that are not in both lists
    for elm in mlist5:  # elm = element
        if elm not in common:
            mseed_list.remove(mseed_dir + elm +
                              '.MSEED')
            print('No picks for %s.' % elm)
    for elem in flist5:  # elem = element
        if elem not in common:
            quakeml_list.remove(xml_dir + elem +
                                '.xml')
            print('No picks for %s.' % elem)

    if len(quakeml_list) != len(mseed_list):
        print('Mseed database does not match Quakeml database!')
        sys.exit()
    elif len(quakeml_list) == len(mseed_list):
        print('Mseed & Quakeml lists match: Continue')

    # loop over ipythoneach QuakeML file and mseed fileevent[0]
    for eventxml, eventmseed in zip(quakeml_list, mseed_list):
        print 'Reading event: %s ' % eventxml
        event = read_events(eventxml)
        # work around, autopick events have no initial origin
        picktime = event[0].picks[0].time
        event[0].origins.append(obspy.core.event.Origin())
        event[0].origins[0].time = picktime
        event[0].origins[0].longitude = 170
        event[0].origins[0].latitude = -40
        event[0].origins[0].depth = 10
        event[0].origins[0].creation_info = ''
        event[0].origins[0].creation_info.version = '_'
        # setting evaluation mode to automatic
        for pick in event[0].picks:
            pick.evaluation_mode = 'automatic'
        # sfile naming
        yr = str(eventxml.split('/')[-1][0:4])
        mo = str(eventxml.split('/')[-1][5:7])
        da = str(eventxml.split('/')[-1][8:10])
        hr = str(eventxml.split('/')[-1][11:13])
        mn = str(eventxml.split('/')[-1][14:16])
        sc = str(eventxml.split('/')[-1][17:19])
        sfilename = da + '-' + hr + mn + '-' + sc + 'L.S' + yr + mo
        write_select(event, sfile_dir + sfilename, 'KPIC', 'L',
                     [eventmseed.split('/')[-1]])
    return


def autolocate_kpick_output(sfile_dir, final_loc_dir, few_picks_dir,
                            P_num, S_num, EQ_RMS):
    """
    Function to autolocate sfiles using SEISAN's hypocenter. The autolocation
    is carried out iteratively removing the picks with the highest residuals
    until some quality criteria are met.

    Arguments:
        sfile_dir {str} -- Directory that contains sfiles for location
        final_loc_dir {str} -- Directory to store the final locations
        few_picks_dir {str} -- Directory to store the locations that failed to
                               meet the criteria
        P_num {int} -- Number of P-wave phases
        S_num {int} -- Number of S-wave phases
        EQ_RMS {float} -- Minimum RMS
    """
    sfile_dir='/Users/home/michaiko/test/sfiles/'
    sfile_list = []
    sfile_list = glob.glob(sfile_dir + '/*')
    sfile_list.sort()

    for sfile in sfile_list:
        print('Locating sfile %s' % sfile.split('/')[-1])
        returned_str = subprocess.check_output(['hyp', sfile])
        if 'no result' in returned_str:
            print(returned_str)
            continue
        try:
            event = read_nordic('hyp.out')
        except Exception as e:
            print(e)
            print('Skipping this sfile')
            continue
        wavefiles = readwavename('hyp.out')

        failed = False
        try:
            while event[0].origins[0].quality.standard_error > EQ_RMS:
                arrival_timeres = [(arrival.pick_id,
                                   abs(arrival.time_residual))
                                   for arrival in event[0].origins[0].arrivals
                                   if arrival.time_residual and
                                   arrival.time_weight == 0]
                arrival_timeres.sort(key=lambda tup: tup[1])
                pick_id_to_remove = arrival_timeres[-1][0]

                # remove the pick
                print('Removing pick...')
                event[0].picks = [pick for pick in event[0].picks
                                  if not pick.resource_id == pick_id_to_remove]

                P_picks_left = [arrival.phase for arrival in
                                event[0].origins[0].arrivals
                                if arrival.time_weight == 0
                                and arrival.phase == 'P']
                S_picks_left = [arrival.phase for arrival in
                                event[0].origins[0].arrivals
                                if arrival.time_weight == 0
                                and arrival.phase == 'S']
                stations = [pick.waveform_id.station_code
                            for pick in event[0].picks]
                stations = list(set(stations))

                # Break the while loop if the location has less than a given
                # number of P-wave picks and/or S-wave picks or picks from
                # different seismic sites.
                if len(P_picks_left) <= P_num and len(S_picks_left) <= S_num:

                    print('Location has {} P-wave picks and {} S-wave picks'
                          ', exiting...'.format(len(P_picks_left),
                                                len(S_picks_left)))
                    failed = True
                    shutil.move('hyp.out', few_picks_dir + '/' +
                                sfile.split('/')[-1])
                    break

                write_select(event, 'hyp.out', 'KPIC', 'L', wavefiles)
                returned_str = subprocess.check_output(['hyp', 'hyp.out'])
                if 'no result' in returned_str:
                    warnings.warn('Cannot locate event, trying as regional')
                    write_select(event, outdir + 'hyp.out', 'KPIC',
                                 'R', wavefiles)
                    returned_str = subprocess.check_output(['hyp', 'hyp.out'])
                    if 'no result' in returned_str:
                        print(returned_str)
                        break
                event = read_nordic('hyp.out')
            if not failed:
                shutil.move('hyp.out', final_loc_dir + '/' +
                            sfile.split('/')[-1])
        except Exception as e:
            print(e)
            continue
    return
