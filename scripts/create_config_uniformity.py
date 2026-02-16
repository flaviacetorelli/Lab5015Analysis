#!/usr/bin/env python
import os, re
#import commands
import math, time
import sys
import argparse
import subprocess

#--------> ex:
#python3 create_config_uniformity.py -t 18 -ml PKU_AA_cfg78 -c config_78.00 -ov 3.0 -r 3992
#python3 create_config_uniformity.py -t 18 -ml DM9001_SM691 -c config_87.00 -ov 0.9 -r 4265 
#python3 create_config_uniformity.py -t 18 -ml DM9003_SM800 -c config_84.00 -ov 0.9 -r 4217
#python3 create_config_uniformity.py -t 18 -ml DM9001_SM691 -c config_87.00 -ov 1.2 -r 4259
# ----
cfgFolder = '/afs/cern.ch/work/f/fcetorel/private/work2/TB_CERN_Sept2025/cfg/TOFHIR2C/uniformity/'
# ----
asicREF = 4
channelMapping0 = { 
0  : [13, 17] , 
1  : [15, 16] , 
2  : [14, 18] , 
3  : [10, 19] , 
4  : [12, 20] , 
5  : [8, 23] , 
6  : [11, 21] , 
7  : [5, 26] , 
8  : [9, 22] , 
9  : [7, 27] , 
10 : [6, 28] , 
11 : [2, 31] , 
12 : [3, 30] , 
13 : [0, 24] , 
14 : [1, 25] , 
15 : [4, 29] , 
}


parser = argparse.ArgumentParser(description='This script creates moduleCharacterization cfg and minEnergy')
parser.add_argument("-ml", "--modulelabel",      required=True,  type=str, help="module label")
parser.add_argument("-r",  "--runs",             required=True,  type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-t",  "--temperature",      required=True,  type=str, help="temperature")
parser.add_argument("-ov", "--Vov",              required=True,  type=str, help="overvoltage")
parser.add_argument("-c",  "--config",           required=True,  type=str, help="config number")
parser.add_argument("-e",  "--extraLabel",       required=False, type=str, help="eg: angle or check or whatever")

args = parser.parse_args()

runs = args.runs

if args.extraLabel:
   label = '%s_Vov%.2f_%s_T%sC' %(args.modulelabel, float(args.Vov),args.extraLabel,  args.temperature)
else:
   label = '%s_Vov%.2f_T%sC' %(args.modulelabel, float(args.Vov) , args.temperature)


#---- write min energy ---

temp_min = '%s/minEnergies_%s_TOFHIR2C.txt'%(cfgFolder,args.modulelabel)
if not (os.path.isfile(temp_min)):
   baseMinEnergy = open('%s/minEnergies_base_TOFHIR2C.txt'%cfgFolder, 'r')
   newMinEnergy  = open('%s/minEnergies_%s_TOFHIR2C.txt'%(cfgFolder,args.modulelabel), 'w')

   command = 'cp %s/minEnergies_base_TOFHIR2C.txt %s/minEnergies_%s_TOFHIR2C.txt'%(cfgFolder, cfgFolder, args.modulelabel)

   os.system(command)

# --- write cfg ---- moduleChar
for b in channelMapping0:
    baseCfg = open('%s/moduleCharacterization_base_TOFHIR2C.cfg'%cfgFolder, 'r')
    if args.extraLabel:
       newCfg = open('%s/moduleCharacterization_%s_refbar%s.cfg'%(cfgFolder,label,b), 'w')
       print ('writing \t moduleCharacterization_%s_refbar%s.cfg'%(label,b))
    else:
       newCfg = open('%s/moduleCharacterization_%s_refbar%s.cfg'%(cfgFolder,label,b ), 'w')
       print ('writing \t moduleCharacterization_%s_refbar%s.cfg'%(label, b))
    
        
    for line in baseCfg:
       #if (line.startswith('Vov') and args.Vov not in line):
          #print ('ERROR: missing ov in moduleCharacterization.cfg file')
          #newCfg.write(line + '%s \n'%args.Vov) # non funziona perche va a capo
          #sys.exit()
       if 'runNumbers' in line:
          newCfg.write(line.replace('runNumbers', '%s'%runs))
       elif 'generalLabel' in line:
          newCfg.write(line.replace('generalLabel', '%s_refbar%s'%(label,b)))
       elif 'moduleLabel' in line:
          newCfg.write(line.replace('moduleLabel', '%s'%args.modulelabel))
       elif 'confNumber' in line:
          newCfg.write(line.replace('confNumber', '%s'%args.config))
          print ('config : ', args.config)
       elif 'vovLabel' in line:
          newCfg.write(line.replace('vovLabel', '%s'%args.Vov))
       elif 'chL' in line:
          newCfg.write(line.replace('CHL', str(asicREF*32+channelMapping0[b][0])))
       elif 'chR' in line: 
          newCfg.write(line.replace('CHR', str(asicREF*32+channelMapping0[b][1])))
       else:
          newCfg.write(line)

baseCfg.close()
newCfg.close()


