
# change the distribution and the geometry
PAR_DIST="EquidistantDistribution"
VOR_DIST="Disk"

# cd ../runs
DIR="${PAR_DIST}_${VOR_DIST}"

mkdir "${DIR}/"

sed -i '' -e"s/VortexInCellManager<T, Dim, .*, .*>/VortexInCellManager<T, Dim, ${PAR_DIST}, ${VOR_DIST}>/" ../alvine/VortexInCell.cpp

cd ../build_serial/alvine/
make

echo "Running ${PAR_DIST} with ${VOR_DIST}..."
./VortexInCell 128 128 500 FFT --overallocate 2.0 --info 10 > ../../runs/$DIR/output.log

# move files
mv particles.csv ../../runs/$DIR/particles.csv
mv energy.csv ../../runs/$DIR/energy.csv

cd ../../plotting
LOCAL_PATH="../runs/${DIR}"

python plot_energy.py $LOCAL_PATH
python plot_initial_distribution.py $LOCAL_PATH
python particle_positions.py $LOCAL_PATH

cd ../runs