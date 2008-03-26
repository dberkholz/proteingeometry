#!/bin/bash

PDBSELECT="/home/donnie/splicer/data/recent.pdb_select90"
RESOLUTION=${1}
PDB_FIELD="2"
RESOLUTION_FIELD="4"
PDB_COUNT="50"
PDB_URL="http://www.pdb.org/pdb/files/"
SCRIPT="get_geom_rmsd.py"
SCRIPT_ARGS="-pe"
# Which lengths/angles do we care about?
GEOM_TYPES=(
    "a3"
    "ome"
    )

# Use local copies, not installed copies
export PYTHONPATH="/home/donnie/src/pgd-utils:${PYTHONPATH}"

if [[ -z ${RESOLUTION} ]]; then
    echo "No resolution argument passed"
    echo "Pass an argument to use as the resolution range."
    echo "Example: 1.0, 2.0, 1.5"
    exit 1
fi

select_pdbs() {
    declare -a PDBS
    local PDBSELECT=$1 RESOLUTION=$2 PDB_COUNT=$3
    local RESOLUTION_FIELD=$4 PDB_FIELD=$5
    # Grab PDB ID and resolution from PDBSelect
    # Select by the resolution we want
    #  Note: be careful with 1.00, because '-1.00' also exists for NMR
    # Sort so we grab the highest-resolution matches
    # Select a certain number of PDBs
    # Convert to lowercase so we can fetch them
    # Drop the chain ID
    PDBS=($(awk '{print $'${RESOLUTION_FIELD}' " " $'${PDB_FIELD}'}' \
        ${PDBSELECT} \
        | grep ^${RESOLUTION} \
        | sort -n \
        | head -n ${PDB_COUNT} \
        | awk '{print $2}' \
        | tr [:upper:] [:lower:] \
        | sed -e "s:^\(.\{4\}\).*:\1:g"))

    # This would be the spot to add automatic resolution increases,
    # if the count is too small

    echo ${#PDBS[@]} ${PDBS[@]}
    return 0
}

read_pdbremarks() {
    local PDBFILE=$1 line i=0
    local -a PDBREMARKS
    while read line; do
        [[ "${line#REMARK}" != "${line}" ]] || continue
        PDBREMARKS[${i}]="${line}"
        (( i++ ))
    done \
        < ${PDBFILE}
    printf "%s\n" "${PDBREMARKS[@]}"
    OLDIFS=$IFS
    IFS='\n'
    echo "${PDBREMARKS[*]}"
    IFS=$OLDIFS
    return 0
}

check_xray() {
    local METHOD
    local PDBFILE=$1
    shift
    for arg in "${@}"; do
        if [[ "${arg#REMARK 200  EXPERIMENT TYPE}" != "${arg}" ]]; then
            METHOD=$(echo "${arg}" | cut -d: -f2)
            break
        fi
    done

    # We only trust X-ray
    if [[ ${METHOD} != *X-RAY\ DIFFRACTION* ]]; then
        echo " Not an X-ray structure: '${METHOD}'" > /dev/stderr
        return 1
    fi
    return 0
}

fetch_pdb() {
    local PDBFILE=$1 URL=$2
    PDBFILEGZ=${PDBFILE}.gz

    if [[ ! -e ${PDBFILE} && ! -e ${PDBFILEGZ} ]]; then
        wget --quiet ${URL}/${PDBFILEGZ}
    fi

    if [[ ! -e ${PDBFILE} && -e ${PDBFILEGZ} ]]; then
        gunzip ${PDBFILEGZ}
    fi

    [[ -e ${PDBFILE} ]] && return 0
    return 1
}

get_value() {
    local VAL=${1##*:}
    VAL=${VAL## }
    VAL=${VAL%% }
    echo ${VAL}
}

get_rfactor() {
    # Keep trying different R factors till we get one
    local RFACTOR PDBFILE=$1
    shift
    for arg in "${@}"; do
        if [[ "${arg}" \
            =~ ^REMARK\ {3}3\ {3}R\ VALUE.*\(WORKING\ SET\) ]]; then
            RFACTOR=$(get_value "${arg}")
            [[ ${RFACTOR} == NULL ]] && continue
            break
        elif [[ "${arg}" \
            =~ ^REMARK\ {3}3\ {3}R\ VALUE.*\(WORKING\ SET,\ NO\ CUTOFF\) ]]; then
            RFACTOR=$(get_value "${arg}")
            [[ ${RFACTOR} == NULL ]] && continue
            break
        elif [[ "${arg}" \
            =~ ^REMARK\ {3}3\ {3}R\ VALUE.*\(WORKING\ \+\ TEST\ SET,\ NO\ CUTOFF\) ]]; then
            RFACTOR=$(get_value "${arg}")
            [[ ${RFACTOR} == NULL ]] && continue
            break
        fi
    done

    echo ${RFACTOR}
    return 0
}

get_freerfactor() {
    # Keep trying different R factors till we get one
    local FREERFACTOR PDBFILE=$1
    shift
    for arg in "${@}"; do
        if [[ "${arg}" \
            =~ ^REMARK\ {3}3\ {3}FREE\ R\ VALUE\ *: ]]; then
            FREERFACTOR=$(get_value "${arg}")
            [[ ${FREERFACTOR} == NULL ]] && continue
            break
        elif [[ "${arg}" \
            =~ ^REMARK\ {3}3\ {3}FREE\ R\ VALUE\ *\(NO\ CUTOFF\) ]]; then
            FREERFACTOR=$(get_value "${arg}")
            [[ ${FREERFACTOR} == NULL ]] && continue
            break
        fi
    done

    echo ${FREERFACTOR}
    return 0
}

get_program() {
    local PROGRAM PDBFILE=$1
    shift
    for arg in "${@}"; do
        if [[ "${arg#REMARK   3   PROGRAM}" != "${arg}" ]]; then
            PROGRAM=$(get_value "${arg}")
            break
        fi
    done
    echo ${PROGRAM}
    return 0
}

get_rmsds() {
    local SCRIPT=$1 SCRIPT_ARGS=$2 GEOM=$3 PDBFILE=$4
    local RESNUM RMSD_EH RMSD_PGD line
    while read -a line; do
        case ${line[0]} in
            "Using")
                RESNUM=${line[1]}
                ;;
            ${GEOM})
                case ${line[4]} in
                    "E&H")
                        RMSD_EH=${line[6]}
                        ;;
                    "PGD")
                        RMSD_PGD=${line[6]}
                        ;;
                    *)
                        echo " Ignoring '${line[4]}'" > /dev/stderr
                        ;;
                esac
                ;;
            '')
                # Skip blank lines
                ;;
            *)
                echo " Ignoring '${line[0]}'" > /dev/stderr
                ;;
        esac
    done \
        < <(${SCRIPT} ${SCRIPT_ARGS} ${GEOM} ${PDBFILE} 2>/dev/null)

    if [[ -z ${RESNUM} || -z ${RMSD_EH} || -z ${RMSD_PGD} ]]; then
        return 1
    fi
    echo ${RESNUM} ${RMSD_EH} ${RMSD_PGD}
    return 0
}

main() {
    local GEOM FILE
    local PDB_ID METHOD PROGRAM RFACTOR RESNUM RMSD_EH RMSD_PGD PDBFILE
    local RMSDS_ARGC STATS
    local -a PDBREM
    set $(select_pdbs \
        ${PDBSELECT} ${RESOLUTION} ${PDB_COUNT} \
        ${RESOLUTION_FIELD} ${PDB_FIELD})
    PDB_COUNT=$1
    shift

    echo "PDB ID count: ${PDB_COUNT}"

    # Clean up from previous runs
    for GEOM in ${GEOM_TYPES[@]}; do
        FILE="${GEOM}_${RESOLUTION}.txt"
        if [[ -e ${FILE} ]]; then
            mv ${FILE} ${FILE}.old
        fi
    done

    # Fetch PDB, get stats
    for PDB_ID in $@; do
        # Initialize to make sure last run's numbers don't stick around
        METHOD=0
        PROGRAM=0
        RFACTOR=0
        RESNUM=0
        RMSD_EH=0
        RMSD_PGD=0

        PDBFILE=${PDB_ID}.pdb
        fetch_pdb ${PDBFILE} ${PDB_URL} || continue

        # Read the PDB into an array, one line per element, for later scanning
        # by check_xray(), get_rfactor(), get_freerfactor(), get_program()
        local arg i=0
        while read arg; do
            PDBREM[$i]="${arg}"
            (( i++ ))
        done < <(read_pdbremarks ${PDBFILE})

        check_xray "${PDBREM[@]}" || continue
        RFACTOR=$(get_rfactor "${PDBREM[@]}")
        FREERFACTOR=$(get_freerfactor "${PDBREM[@]}")
        PROGRAM=$(get_program "${PDBREM[@]}")

        # Run script, save values (to file?)
        for GEOM in ${GEOM_TYPES[@]}; do
            RMSDS_ARGC=${#GEOM_TYPES[@]}
            # Add one for RESNUM
            ((RMSDS_ARGC++))

            set $(get_rmsds ${SCRIPT} ${SCRIPT_ARGS} ${GEOM} ${PDBFILE}) \
                || continue

            # Ignore this PDB if we couldn't calculate stats for it
            if (( $# != RMSDS_ARGC )); then
                echo " Ignoring ${PDB_ID}" > /dev/stderr
                break
            fi
            RESNUM=$1
            RMSD_EH=$2
            RMSD_PGD=$3
            STATS="\"${PDB_ID}\"\t${RESNUM}\t\"${PROGRAM}\"\t${RFACTOR}\t${FREERFACTOR}\t${RMSD_EH}\t${RMSD_PGD}"
            echo -e "${STATS}"
            echo -e "${STATS}" >> ${GEOM}_${RESOLUTION}.txt
        done
    done
}

main
