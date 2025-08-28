#!/usr/bin/env bash
#* This file is part of the MOOSE framework
#* https://mooseframework.inl.gov
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html
set -ex
echo "[DEBUG] USING AIKUBO FIX"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -n "$WASP_SRC_DIR" ]; then
  SKIP_SUBMODULE_UPDATE=1
else
  MOOSE_DIR=${SCRIPT_DIR}/..
  WASP_SRC_DIR=${MOOSE_DIR}/framework/contrib/wasp
fi

# Loop over specified command line arguments and turn on any found options
for ARG in "$@" ; do
  if   [[ "${ARG}" == "--help" ]]                  ; then HELP=1
  elif [[ "${ARG}" == "--fast" ]]                  ; then FAST=1
  elif [[ "${ARG}" == "--skip-submodule-update" ]] ; then SKIP_SUBMODULE_UPDATE=1
  else                                               EXTRA_ARGS+=("$ARG")
  fi
done

# Display the help menu to show the list of available options if requested
if [[ -n "${HELP}" ]] ; then
  echo "Usage: $0 [ --help | --fast | --skip-submodule-update ]"
  echo
  echo "--help                   Display this message listing the options that are available"
  echo "--skip-submodule-update  Do not update the WASP submodule, use the current version"
  echo "--fast                   Run WASP 'make install' only and do NOT update or configure"
  echo "*********************************************************************"
  echo
  exit 0
fi

# If we are going fast then check and make sure the build directory exists
if [[ -n "${FAST}"  ]] ; then
  if [[ ! -d "${WASP_SRC_DIR}"/build ]] ; then
    echo "Error: a build directory must exist to use the --fast option"
    exit 1
  fi
  cd "${WASP_SRC_DIR}"/build
# If we are not going fast then update WASP, remove build, and reconfigure
else
  if [ -z "$SKIP_SUBMODULE_UPDATE" ]; then
    cd "$MOOSE_DIR"
    echo "[DEBUG]"
    echo "tree -a -L 3"
    git_dir=$(git rev-parse --show-cdup 2>/dev/null)
    if [[ $? -eq 0 ]] && [ -z "$SKIP" ] && [[ "$git_dir" == "" ]]; then
      git submodule update --init --recursive "${WASP_SRC_DIR}"
    fi
  fi

  rm -rf "${WASP_SRC_DIR}"/build
  mkdir -p "${WASP_SRC_DIR}"/build
  cd "${WASP_SRC_DIR}"/build
  WASP_OPTIONS="-DCMAKE_INSTALL_PREFIX:STRING=${WASP_PREFIX:-${WASP_SRC_DIR}/install}"
  source $SCRIPT_DIR/configure_wasp.sh
  configure_wasp "$WASP_OPTIONS" ../ "${EXTRA_ARGS[@]}"
  if [[ $? -ne 0 ]] ; then
    echo "Error: configure step for WASP failed to complete successfully"
    exit 1
  fi
fi

# Build with MOOSE_JOBS jobs if defined, otherwise build with a single job
make -j ${MOOSE_JOBS:-1} install
if [[ $? -ne 0 ]] ; then
  echo "Error: build step for WASP failed to complete successfully"
  exit 1
fi
echo "[DEBUG] wasp configure successful"
exit 0
