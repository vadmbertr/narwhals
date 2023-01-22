#!/bin/bash

# Remplacer submission_example par la vraie version (qui correspond au nom d'une branche)
# Enlever les # qui trainent à la fin des lignes commençant par #OAR

#OAR -n submission_example                    # nom du job, A CHANGER
#OAR -l /nodes=1/core=10,walltime=01:00:00    # ressources, A CHANGER
#OAR --stdout submission_example_log.out      # fichier de sortie, A CHANGER
#OAR --stderr submission_example_log.err      # fichier de sortie (erreurs), A CHANGER
#OAR --project pr-whales                      # nom du projet

version="submission_example"                  # version / branche, A CHANGER

source /applis/site/nix.sh

folder_name="narwhals_${version}"   # dossier où est cloné la branche du repo git
if [ -d "${folder_name}" ] 
then
    rm -rf "${folder_name}"         # supprimé si existe déjà
fi
mkdir "${folder_name}"              # création

save_folder="/bettik/PROJECTS/pr-whales/COMMON/"${version}  # dossier pour sauvegarder les sorties du script R
mkdir -p "${save_folder}"

cd "${folder_name}"
git clone -b "${version}" --single-branch git@github.com:vadmbertr/narwhals.git .
cd R
# lance le script, CHANGER les arguments
Rscript 1_"${version}".R /bettik/PROJECTS/pr-whales/COMMON/db_narwhal_2017_2018.txt arg2 arg3 "${save_folder}"
