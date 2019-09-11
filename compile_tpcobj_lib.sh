#!/bin/bash

# This script will compile the libraries in Modules
# and move them to the lib area
BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MODULETOP=${PWD}
MODULEANA=${MODULETOP}/Modules/

# Now we can build the xsecAna classes and create the libraries and copy to lib area
echo -e "${BLUE}Compiling TPC Object libraries...${NC}"
cd ${MODULEANA}
make clean
make
make copy

cd ${MODULETOP}
echo -e "${BLUE}You are now in: " ${MODULETOP}
echo -e "${GREEN}Setup Finished ${NC}"
