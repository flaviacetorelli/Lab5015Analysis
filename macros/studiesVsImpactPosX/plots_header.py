#! /usr/bin/env python
import math

import ROOT

def draw_logo():
    logo_x = 0.17
    logo_y = 0.94
    logo = ROOT.TLatex()
    logo.SetNDC()
    logo.SetTextSize(0.045)
    logo.SetTextFont(62)
    logo.DrawText(logo_x , logo_y,'CMS')
    logo.SetTextFont(52)
    logo.DrawText(logo_x + 0.08, logo_y, 'Preliminary')
    logo.SetTextFont(62)
    logo.DrawText(0.80,0.94,'Phase-2')
    
    return logo
