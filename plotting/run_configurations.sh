
# change the distribution and the geometry
PAR_DIST="EquidistantDistribution"
VOR_DIST="JetPenetration"
TIME_STEPS=1000
GRID_SIZE=128
SOLVER="FFT"

cd ../runs
DIR="${PAR_DIST}_${VOR_DIST}"
LOCAL_PATH="../runs/${DIR}"

mkdir "${DIR}/"

sed -i '' -e"s/VortexInCellManager<T, Dim, .*, .*>/VortexInCellManager<T, Dim, ${PAR_DIST}, ${VOR_DIST}>/" ../alvine/VortexInCell.cpp

cd ../build_serial/alvine/
make

echo "Running ${PAR_DIST} with ${VOR_DIST}..."
./VortexInCell $GRID_SIZE $GRID_SIZE $TIME_STEPS $SOLVER --info 10 > ../../runs/$DIR/output.log

# move files
mv particles.csv ../../runs/$DIR/particles.csv
mv energy.csv ../../runs/$DIR/energy.csv

cd ../../plotting

python plot_energy.py $LOCAL_PATH
python plot_initial_distribution.py $LOCAL_PATH
python plot_4_frames.py $LOCAL_PATH
python particle_positions.py $LOCAL_PATH