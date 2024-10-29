echo "begin very strange procedure - mc16d"

#run_ntuple --config bkg_myy90_175_dr0p5.config &
#run_ntuple --config bkg_myy90_175_dr1p0.config &

run_ntuple --config gamgam.config &
#run_ntuple --config bkg_Wplus_dr1p0.config &
#run_ntuple --config bkg_Wplus_dr1p5.config &

wait

echo "eeeend!!"
