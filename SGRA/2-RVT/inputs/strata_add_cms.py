import json
import numpy as np
import copy

def strata_add_cms(input_fn, cms_p, cms_sa, damping=5, magnitude=None, distance=None, duration=None, description=None, name=None, out_fn = None):
    with open(input_fn, 'r') as f:
        data = json.load(f)
        motion_library = data['motionLibrary']
        motions = motion_library['motions']
        cur_motion = motions[-1]
        # print('cur_motion pga=', cur_motion['targetRespSpec']['sa'][0])
        new_motion = copy.deepcopy(cur_motion)
        
        # new_motion = cur_motion.copy()
        # print('new_motion pga=', new_motion['targetRespSpec']['sa'][0])
        if not magnitude is None:
            new_motion['magnitude'] = magnitude
        if not duration is None:
            new_motion['duration'] = duration
        if not distance is None:
            new_motion['distance'] = distance
        if not description is None:
            new_motion['description'] = description
        if not name is None:
            new_motion['name'] = name
        if out_fn is None:
            out_fn = input_fn
        new_motion['targetRespSpec']['damping'] = damping
        new_motion['targetRespSpec']['period'] = list(cms_p)
        new_motion['targetRespSpec']['sa'] = list(cms_sa)
        motions.append(new_motion)
        # print('cur_motion pga2=', cur_motion['targetRespSpec']['sa'][0])
        # print('new_motion pga2=', new_motion['targetRespSpec']['sa'][0])
    with open(out_fn, 'w') as f:
        json.dump(data,f, indent=4)
    # print('cur_motion pga3=', cur_motion['targetRespSpec']['sa'][0])
    # print('new_motion pga3=', new_motion['targetRespSpec']['sa'][0])
    # exit()

if __name__ == '__main__':
    
    ##input json file
    json_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\SGRA\2-RVT\inputs/'
    json_file = 'Salem_RVT_empirical_inputs.json'
    ## read deagg results
    deagg_file = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\SGRA\0-deaggResults/Salem_deagg_summary.csv'
    deaggdata  = np.genfromtxt(deagg_file,delimiter=',',dtype=None, names=True)
    # x = cmsdata
    deaggdata = deaggdata.view((float, len(deaggdata.dtype.names)))
    
    # read the CMS periods and spectral accelerations
    cms_folder = r'C:\Users\FEli\Golder Associates\20368492, National Grid Lynn LNG MA - 5 Technical Work\SGRA\1-CMS/'
    list_rp = [50, 100, 250, 475, 975, 2475, 5000, 10000]
    # list_rp = [50,100]
    list_target = ['uhs', 'uhs', 'cms','cms','cms',
                   'cms','cms','cms']
    list_periods = [[0.2], [0.2], [0.2,1], [0.2,1,5], [0.2,1,5], 
                    [0.2,1,5], [0.2,1,5], [0.2,1,5], ]
    cms_initial = [[0.01,0.02,0.03,0.05,0.075,
                0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,
                1,1.5,2,3,4,5,7.5,10]]
    cases = 'perids_s'
    for ii in range(len(list_rp)):
        irp = list_rp[ii]
        cms_file = cms_folder + f'CMS_Salem_{irp}yr.csv'
        cmsdata  = np.genfromtxt(cms_file, delimiter=',', dtype=None, names=True)
        x = cmsdata
        x = x.view((float, len(x.dtype.names)))
        
        periods = list_periods[ii]
        target = list_target[ii]
        
        header = np.char.strip(cmsdata.dtype.names)
   
        # print('deagg', deaggdata[:,0])
        cmsdata = x
        cms_p = cmsdata[:,0].tolist()
        
        for ip in periods:
            temp = np.where((deaggdata[:,0]==irp)&(deaggdata[:,1]==ip))
            # print(temp, ' ip=', ip, ' irp=',irp)
            # print(deaggdata[temp][0])
            
            mag = round(deaggdata[temp][0][2],2)
            distance = round(deaggdata[temp][0][3],2)
            d75 = round(deaggdata[temp][0][4],2)
            if target == 'uhs':
                case = "target"
            else:
                # mag_int = int(mag*100)
                # distance_int = int(distance*100)
                case = f'Mw{mag}_R{distance}_T{ip}' #  "Mw5.84_R174.6_T1"
                case = case.replace('.','')
                
            index = np.where(header==case)
            print('###case', case, ' mag=', mag, ' dist=', distance)
            # print('header=', header, ', index=', index)
            
            cms_sa = list(cmsdata[:,index].flat)
            print('pga = ', cms_sa[0])
            cms_initial.append(cms_sa)
            cases += f',yr{irp}+p{ip}s'
            # print('type(cms_sa', type(cms_sa))
            # print(cms_sa)
            # print("###")
            # print(list(cms_sa))
            # print("index=", index)
            # print(type(cms_p))
            # cms_p = np.logspace(-3,1,50)
            # cms_sa = 0*cms_p + 1
            # make sure the json file does not include the Strata results
            strata_add_cms(json_folder+json_file, cms_p, cms_sa, 
                           name=f'Comptible RVT Motion (M{mag} @ {distance}km)', 
                           damping=5, magnitude=mag, distance=distance, duration=d75, 
                           description=f'{irp}yr at {ip}s', out_fn=json_file)
    
    cms_initial = np.transpose(cms_initial)
    fname = json_folder+ 'Salem_cms_initial.csv' #'CMS_NGA_east.csv'
    with open(fname, 'w') as f:    
        # np.savetxt(f, case, delimiter=',')
        f.write(cases + '\n')
    with open(fname, 'ba') as f:
        np.savetxt(f, cms_initial, delimiter=',')
    # return()