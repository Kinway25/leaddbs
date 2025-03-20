'''
    By K. Butenko
    Runs OSS-DBS to compute 'VTRs' for sEEG
'''

import pandas as pd
import numpy as np
import sys
import os
import json
import subprocess
import re

from run_OSS4SEEG_Stim_no_shift import check_electrode_availability, get_geom_definitions
from ossdbs.electrodes.defaults import default_electrode_parameters

if __name__ == '__main__':

    # called from MATLAB
    # sys.argv[1] - full path to the reconstruction file
    # sys.argv[2] - 'CC' if current-controlled
    # sys.argv[3] - Electrode ID (as integer), optional 

    ''' input processing '''
    SEEG_recos = sys.argv[1]
    _,extension = os.path.splitext(SEEG_recos)
    if extension == '.tsv':
        # either we get reco from tsv (Clemens' format)
        SEEG_recos_df = pd.read_csv(SEEG_recos, sep='\t')
    elif extension == '.mat':
        print("Lead-DBS reconstruction files are currently not supported")
        raise SystemExit()
        # reconstruction (a la Garance)
        import h5py
        SEEG_recos_mat = h5py.File(str(SEEG_recos), "r")
        # create a pandas dataframe analogous to above

    if sys.argv[2] == 'CC':
        current_controlled = True
    else:
        current_controlled = False

    if len(sys.argv) > 3:
        Electrode_ID = int(sys.argv[3])
        SEEG_recos_df = SEEG_recos_df[SEEG_recos_df['Electrode_ID'] == Electrode_ID]
    else:
        Electrode_ID = None

    # some auto-definitions
    stim_folder = os.path.dirname(SEEG_recos)
    Amplitude = [0.001]  # 1 mA. The exact value is not important for VTRs
    contacts2simulate = SEEG_recos_df.name   # all reconstructed

    # we iterate over all contacts of all electrodes, but rebuilding the geometry everytime
    for cnt_i in range(len(contacts2simulate)):
    
        # electrode type and ID
        oss_electrode = check_electrode_availability(SEEG_recos_df['electrode'][cnt_i])  
        if Electrode_ID == None:
            # definition from the reconstruction sheet
            Electrode_ID = SEEG_recos_df['Electrode_ID'][cnt_i]
    
        ''' determine two contacts (active and adjacent) to build the trajctory '''
        cnt_ID = contacts2simulate[cnt_i]  # actual label
        flip = False;
        
        # IMPORTANT: this is a hard assumption that contact labels start with 1!
        index_on_electrode = int(re.sub(r"\D", "",cnt_ID)) - 1  # integer index
        if index_on_electrode < 0:
            print("Numbering for contact labels is expected to start from 1!")
            raise SystemExit
            
        if cnt_i+1 == len(contacts2simulate):
            # last listed contact
            cnt_ID2 = contacts2simulate[cnt_i-1]
            flip = True
        else:
            cnt_ID2 = contacts2simulate[cnt_i+1]
            index_on_electrode2 = int(re.sub(r"\D", "",cnt_ID2)) - 1
            if index_on_electrode2 - index_on_electrode != 1:
                # last contact, use previous to define a trajectory
                cnt_ID2 = contacts2simulate[cnt_i-1]
                flip = True
        
        # find indices of these contacts in SEEG_recos_df
        inx = SEEG_recos_df.index[SEEG_recos_df['name'] == cnt_ID]
        # check if the second contact coordinates are present
        if (SEEG_recos_df['name'].eq(cnt_ID2)).any():
            inx2 = SEEG_recos_df.index[SEEG_recos_df['name'] == cnt_ID2]
        else:
            print("The second contact was not found, cannot reconstruct the trajectory!")         
            raise SystemExit
        
        # these are two contact coordinates to define the trajectory
        contact_coords = [[SEEG_recos_df['x'][inx].values[0], SEEG_recos_df['y'][inx].values[0], SEEG_recos_df['z'][inx].values[0]],[SEEG_recos_df['x'][inx2].values[0], SEEG_recos_df['y'][inx2].values[0], SEEG_recos_df['z'][inx2].values[0]]]
        
        '''  now get the implantation coordinates '''
        Dimensions,grid_center,unit_directions = get_geom_definitions(contact_coords)
        if flip:
            # flip if the second contact is distal
            unit_directions = unit_directions * -1.0
 
        elec_params = default_electrode_parameters[oss_electrode]
        # this might be wrong if the first contact (active tip) has a different length
        imp_coords = np.array([contact_coords[0][0],contact_coords[0][1],contact_coords[0][2]]) - index_on_electrode * (elec_params.contact_length + elec_params.contact_spacing) * unit_directions  
        
        # offset = from tip to the center of the first contact 
        offset = elec_params.get_center_first_contact() * 1.0
        tip_position = imp_coords - offset * unit_directions
            
        for amp in Amplitude:
        
            custom_params = {
    
                "BrainRegion": {
                    # center at the head marker
                    "Center": {
                        "x[mm]":  grid_center[0],
                        "y[mm]":  grid_center[1],
                        "z[mm]":  grid_center[2],
                    },
                    "Dimension": {
                      "x[mm]": Dimensions['x[mm]'],
                      "y[mm]": Dimensions['y[mm]'],
                      "z[mm]": Dimensions['z[mm]']
                    },
                    "Shape": "Ellipsoid"
                },
                "Electrodes": [
                    {
                      "Name": oss_electrode,
                      "Rotation[Degrees]": -179.18978179152586,  # not relevant
                      "Direction": {
                        "x[mm]": unit_directions[0],
                        "y[mm]": unit_directions[1],
                        "z[mm]": unit_directions[2]
                      },
                      "TipPosition": {
                        "x[mm]": tip_position[0],
                        "y[mm]": tip_position[1],
                        "z[mm]": tip_position[2]
                      },
                      "Contacts": [         # we need only one contact! Move the forth one for now!
                        {
                          "Contact_ID": index_on_electrode+1,
                          "Active": False,
                          "Current[A]": amp,
                          "Voltage[V]": False,
                          "Floating": True,
                          "SurfaceImpedance[Ohmm]": {
                            "real": 0.0,
                            "imag": 0.0
                          },
                          "MaxMeshSize": 1000000.0,
                          "MaxMeshSizeEdge": 0.08545132017764238  # hardcoded for now
                        }
                      ],                      
                      "EncapsulationLayer": {
                        "Thickness[mm]": 0.0,
                        "Material": "White matter",
                        "DielectricModel": "ColeCole4",
                        "DielectricParameters": None,
                        "MaxMeshSize": 1000000.0
                      }
                    }
                  ],
    
                "Surfaces": [
                    {
                        "Name": "BrainSurface",
                        "Active": True,
                        "Current[A]": -1*amp,
                        "Voltage[V]": 0.0,
                    }
                ],
                "MaterialDistribution": {
                    "MRIPath": os.path.join(stim_folder,'segmask.nii'),
                    "DiffusionTensorActive": False,
                    "DTIPath": '',
                },
                "DielectricModel": {
                    "Type": "ColeCole4",
                    "CustomParameters": None
                },
                "Mesh": {
                  "LoadMesh": False,
                  "LoadPath": "",
                  "MeshingHypothesis": {
                    "Type": "Default",
                    "MaxMeshSize": 1000000.0,
                    "MeshSizeFilename": ""
                  },
                  "HPRefinement": {
                    "Active": False,
                    "Levels": 2,
                    "Factor": 0.125
                  },
                  "AdaptiveMeshRefinement": {
                    "Active": False,
                    "MaxIterations": 10,
                    "ErrorTolerance": 0.1
                  },
                  "MaterialRefinementSteps": 1,
                  "MeshSize": {
                    "Edges": {},
                    "Faces": {},
                    "Volumes": {}
                  },
                  "SaveMesh": False,
                  "SavePath": "mesh"
                },
                "EQSMode": False,
                "FEMOrder": 2,
                "ComputeImpedance": False,
                "StimulationSignal": {
                    "Type": "Multisine",
                    "ListOfFrequencies": [
                        10000
                    ],
                    "Frequency[Hz]": 130.0,
                    "PulseWidth[us]": 60.0,   # not relevant for stim volumes
                    "PulseTopWidth[us]": 0.0,
                    "CounterPulseWidth[us]": 0.0,
                    "InterPulseWidth[us]": 0.0,
                    "SpectrumMode": "FullSpectrum",
                    "CutoffFrequency": 100000000.0,
                    "CurrentControlled": current_controlled
                },
                "Solver": {
                  "Type": "CG",
                  "Preconditioner": "local",
                  "PreconditionerKwargs": {},
                  "MaximumSteps": 1000,
                  "Precision": 1e-10
                },
                "PointModel": {
                    "Lattice": {
                        "Active": True,
                        "Center": {
                            "x[mm]": grid_center[0],
                            "y[mm]": grid_center[1],
                            "z[mm]": grid_center[2]
                        },
                        "Shape": {"x": 51, "y": 51, "z": 51},
                        "Direction": {
                            "x[mm]": 0,
                            "y[mm]": 0,
                            "z[mm]": 1,
                        },
                        "PointDistance[mm]": 0.4,
                        "CollapseVTA": True,  # questionable
                    }
                },
                "OutputPath": os.path.join(os.path.dirname(SEEG_recos),'Results_VTR_E' + str(Electrode_ID) + '_1mA_' + cnt_ID),
                "SaveImpedance": False,
                "ExportVTK": True,
                "TemplateSpace": False,
                "ModelSide": 0,
                "CalcAxonActivation": False,
                "ActivationThresholdVTA[V-per-m]": 200.0,
                "FailFlag": 'rh',
                "OutOfCore": False
            }
    
            # save the settings
            json_settings = json.dumps(custom_params, indent=2)
            with open(os.path.join(stim_folder,'oss-dbs_parameters.json'), "w") as outfile:
                outfile.write(json_settings)
    
            # run OSS-DBS
            with open(os.devnull, 'w') as FNULL: subprocess.call('ossdbs ' + os.path.join(stim_folder,'oss-dbs_parameters.json'),shell=True)
            