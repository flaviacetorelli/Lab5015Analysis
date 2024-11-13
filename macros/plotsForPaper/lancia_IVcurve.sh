# DCR vs Vov Plots for JINST
python3 plot_IVcurve.py --IVorDCR DCR --comparison cellsize    --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsOV_minorRevJINST_Nov24/ #paper
python3 plot_IVcurve.py --IVorDCR DCR --comparison irradiation --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsOV_minorRevJINST_Nov24/ #paper
python3 plot_IVcurve.py --IVorDCR DCR --comparison temperature --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsOV_minorRevJINST_Nov24/
python3 plot_IVcurve.py --IVorDCR DCR --comparison vendor      --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsOV_minorRevJINST_Nov24/
python3 plot_IVcurve.py --IVorDCR DCR --comparison jinst      --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsOV_minorRevJINST_Nov24/ #paper

# DCR vs Gain Plots
python3 plot_DCRvsGain.py --comparison cellsize    --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsGain_minorRevJINST_Nov24/ #paper
python3 plot_DCRvsGain.py --comparison irradiation --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsGain_minorRevJINST_Nov24/
python3 plot_DCRvsGain.py --comparison temperature --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsGain_minorRevJINST_Nov24/
python3 plot_DCRvsGain.py --comparison vendor      --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/DCRvsGain_minorRevJINST_Nov24/
