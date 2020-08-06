rm -rf 0
cp -r 0_org 0

foamCleanPolyMesh

blockMesh

snappyHexMesh -overwrite > snappyHexMesh.log

createPatch -overwrite > createPatch.log





