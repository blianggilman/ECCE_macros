#!/bin/bash

outdir="/project/projectdirs/m3763/blianggi"

if [[ ! -e $outdir/out_ECCE ]]; then
	mkdir $outdir/out_ECCE
fi

if [[ ! -e $outdir/out_ECCE/$1 ]]; then
	mkdir $outdir/out_ECCE/$1
fi

core_dir="/project/projectdirs/m3763/blianggi/ECCE_setup"
source $core_dir/Singularity/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh -n
export MYINSTALL=$core_dir/Singularity/install
source $core_dir/Singularity/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/setup_local.sh $MYINSTALL

export g4_dir=$core_dir/Singularity/cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core
export G4LEVELGAMMADATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/PhotonEvaporation5.5
export G4_MAIN=$g4_dir/geant4.10.06.p02
export G4LEDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4EMLOW7.9.1
export G4NEUTRONHPDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4NDL4.6
export G4ENSDFSTATEDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4ENSDFSTATE2.2
export G4RADIOACTIVEDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/RadioactiveDecay5.4
export G4ABLADATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4ABLA3.1
export G4PIIDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4PII1.3
export G4PARTICLEXSDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4PARTICLEXS2.1
export G4SAIDXSDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4SAIDDATA2.0
export G4REALSURFACEDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/RealSurface2.1.1
export G4INCLDATA=$g4_dir/geant4.10.06.p02/share/Geant4-10.6.2/data/G4INCL1.0

cd $core_dir/ECCE_macros/detectors/EICDetector/
root -b -q "Fun4All_G4_EICDetector.C($3,\"$outdir/out_ECCE/$1\",\"$outdir/out_ECCE/$1/G4EICDetector_$2\")"