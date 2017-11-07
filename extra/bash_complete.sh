_Ant() 
{
    local cur prev opts cmd
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    cmd=$(which ${COMP_WORDS[0]})
    cmddir=$(dirname ${cmd})
    cmdbase=$(basename ${cmd})
    Ant_completion="${cmddir}/Ant-completion"
    if ! [ -f ${Ant_completion} ]; then
      	return 1
    fi 

    case "${prev}" in
        -i|-o|--input)
            _filedir '@(root|xz)'
            return 0
            ;;

        -s|--setup)
            local setups=$( ${Ant_completion} ${cmdbase} --setup)
            COMPREPLY=( $(compgen -W "${setups}" -- ${cur}) )
            return 0
            ;;
        -p|--physics)
            local physics=$( ${Ant_completion} ${cmdbase} --physics )
            COMPREPLY=( $(compgen -W "${physics}" -- ${cur}) )
            return 0
            ;;
        -c|--calibration)
            for (( i = 0 ; i < ${#COMP_WORDS[@]}-1 ; i++ )); do
                word=${COMP_WORDS[${i}]} 
                case "${word}" in
                    -s|--setup)
                        local cals=$( ${Ant_completion} ${cmdbase} --calibration ${COMP_WORDS[$((i+1))]} )
                        COMPREPLY=( $(compgen -W "${cals}" -- ${cur}) )
                        return 0
                        ;;
                esac
            done
            return 0
            ;;
    esac

    opts="--help --version --batch --u_writeuncalibrated --u_disablereconstruct --u_writecalibrated --p_disableParticleID -i --input -s --setup -p --physics -o --output -v --verbose -m --maxevents -O -c --calibration"
    if [[ ${cur} == * ]] ; then
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
        return 0
    fi
}
complete -F _Ant Ant
