_Ant() 
{
    local cur prev opts cmd
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    cmd="${COMP_WORDS[0]}"

    case "${prev}" in
        -i|-o|--input)
            _filedir '@(root|xz)'
            return 0
            ;;

        -s|--setup)
            local setups=$( "${cmd}" --list-setups )
            COMPREPLY=( $(compgen -W "${setups}" -- ${cur}) )
            return 0
            ;;
        -p|--physics)
            local physics=$( "${cmd}" --list-physics )
            COMPREPLY=( $(compgen -W "${physics}" -- ${cur}) )
            return 0
            ;;
        -c)
            local cals=$( "${cmd}" --list-calibrations )
            COMPREPLY=( $(compgen -W "${cals}" -- ${cur}) )
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
