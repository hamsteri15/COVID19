rm -rf 0
cp -r 0_org 0

mapFields -sourceTime "latestTime" ../create_turbulence

buoyantPimpleFoam > buoyantPimpleFoam.log
