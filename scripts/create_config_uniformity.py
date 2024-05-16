#!/usr/bin/env python
import os, re
#import commands
import math, time
import sys
import argparse
import subprocess

#--------> ex: 
#python3 create_config_uniformity.py -d /afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/ -r 6062 -t -35 -ov 1.5 -ml HPK_2E14_C25_LYSO100056 -c config_51.04
#python3 create_config_uniformity.py -d /afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/ -r 6260,6264 -t 5 -ov 3.5 -ml HPK_nonIrr_C25_LYSO818  -c config_61.00
# ----
#cfgFolder = '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/scripts/cfg/TOFHIR2C/'
# ----


channelMapping0 = { 
0  : [14, 17] , 
1  : [12, 19] , 
2  : [10, 21] , 
3  : [8, 23] , 
4  : [6, 25] , 
5  : [4, 27] , 
6  : [2, 29] , 
7  : [0, 31] , 
8  : [1, 30] , 
9  : [3, 28] , 
10 : [5, 26] , 
11 : [7, 24] , 
12 : [9, 22] , 
13 : [11, 20] , 
14 : [13, 18] , 
15 : [15, 16] , 
}


parser = argparse.ArgumentParser(description='This script creates moduleCharacterization cfg and minEnergy')
parser.add_argument("-d",  "--directory",        required=True,  type=str, help="folder of cfgs")
parser.add_argument("-ml", "--modulelabel",      required=True,  type=str, help="module label")
parser.add_argument("-r",  "--runs",             required=True,  type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-t",  "--temperature",      required=True,  type=str, help="temperature")
parser.add_argument("-ov", "--Vov",              required=True,  type=str, help="overvoltage")
parser.add_argument("-c",  "--config",           required=True,  type=str, help="config number")
parser.add_argument("-e",  "--extraLabel",       required=False, type=str, help="eg: angle or check or whatever")
parser.add_argument("-u",  "--uniformity",       action='store_true', help="cfg files to study uniformity")

args = parser.parse_args()



runs = args.runs
cfgFolder = args.directory

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


if args.uniformity: 
    
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
              newCfg.write(line.replace('CHL', str(channelMapping0[b][0])))
           elif 'chR' in line: 
              newCfg.write(line.replace('CHR', str(channelMapping0[b][1])))
           else:
              newCfg.write(line)

    baseCfg.close()
    newCfg.close()

else:

    baseCfg = open('%s/moduleCharacterization_base_TOFHIR2C.cfg'%cfgFolder, 'r')
    if args.extraLabel:
       newCfg = open('%s/moduleCharacterization_%s.cfg'%(cfgFolder,label), 'w')
       print ('writing \t moduleCharacterization_%s.cfg'%(label))
    else:
       newCfg = open('%s/moduleCharacterization_%s.cfg'%(cfgFolder,label), 'w')
       print ('writing \t moduleCharacterization_%s.cfg'%(label))
    
        
    for line in baseCfg:
       if (line.startswith('Vov') and args.Vov not in line):
          print ('ERROR: missing ov in moduleCharacterization.cfg file')
          newCfg.write(line + '%s \n'%args.Vov) # non funziona perche va a capo
          sys.exit()
       elif 'runNumbers' in line:
          newCfg.write(line.replace('runNumbers', '%s'%runs))
       elif 'generalLabel' in line:
          newCfg.write(line.replace('generalLabel', '%s'%label))
       elif 'moduleLabel' in line:
          newCfg.write(line.replace('moduleLabel', '%s'%args.modulelabel))
       elif 'confNumber' in line:
          newCfg.write(line.replace('confNumber', '%s'%args.config))
          print ('config : ', args.config)
       elif 'chL' in line:
          newCfg.write(line.replace('CHL', '1'))
       elif 'chR' in line: 
          newCfg.write(line.replace('CHR', '30'))
       else:
          newCfg.write(line)

    newCfg.close()
    baseCfg.close()
    
#    if args.uniformity:
#        for b in channelMapping0:
#           with open("%s/moduleCharacterization_%s.cfg"%(cfgFolder, label), 'r') as fi:
#               contents = fi.read()
#               replaced_contents = contents.replace("CHL", str(channelMapping0[b][0])).replace("CHR",  str(channelMapping0[b][1]))
#               with open("%s/moduleCharacterization_%s_refbar%i.cfg"%(cfgFolder, label, b), "w") as fo:
#                   fo.write(replaced_contents)

# --- write cfg ---- pulseShape
#baseCfg = open('%s/drawPulseShapeTB_base.cfg'%cfgFolder, 'r')
#
#
#if args.extraLabel:
#   newCfg = open('%s/drawPulseShapeTB_%s.cfg'%(cfgFolder,label), 'w')
#   print 'writing \t drawPulseShapeTB_%s.cfg'%(label)
#else:
#   newCfg = open('%s/drawPulseShapeTB_%s.cfg'%(cfgFolder,label), 'w')
#   print 'writing \t drawPulseShapeTB_%s.cfg'%(label)
#
#
#for line in baseCfg:
#   if (line.startswith('Vov') and args.Vov not in line):
#      newCfg.write(line + '%s'%args.Vov)
#   elif 'runNumbers' in line:
#      newCfg.write(line.replace('runNumbers', '%s'%runs))
#   elif 'generalLabel' in line:
#      newCfg.write(line.replace('generalLabel', '%s'%label))
#   elif 'moduleLabel' in line:
#      newCfg.write(line.replace('moduleLabel', '%s'%args.modulelabel))
#   elif 'confNumber' in line:
#      newCfg.write(line.replace('confNumber', '%s'%args.config))
#      print 'config : ', args.config
#
#   else:
#      newCfg.write(line)
#
#baseCfg.close()
#newCfg.close()




