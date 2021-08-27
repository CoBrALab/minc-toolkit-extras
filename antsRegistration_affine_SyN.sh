#!/bin/bash
# ARG_HELP([The general script's help msg])
# ARG_POSITIONAL_SINGLE([movingfile],[The moving image])
# ARG_POSITIONAL_SINGLE([fixedfile],[The fixed image])
# ARG_POSITIONAL_SINGLE([outputbasename],[The basename for the output transforms])
# ARG_OPTIONAL_SINGLE([moving-mask],[],[Mask for moving image],[NOMASK])
# ARG_OPTIONAL_SINGLE([fixed-mask],[],[Mask for fixed image],[NOMASK])
# ARG_OPTIONAL_SINGLE([resampled-output],[o],[Output resampled file])
# ARG_OPTIONAL_SINGLE([resampled-linear-output],[],[Output resampled file with only linear transform])
# ARG_OPTIONAL_SINGLE([initial-transform],[],[Initial linear transform],[NONE])
# ARG_OPTIONAL_SINGLE([linear-type],[],[Type of affine transform],[affine])
# ARG_OPTIONAL_BOOLEAN([close],[],[Images are starting off close, skip large scale pyramid search],[])
# ARG_OPTIONAL_REPEATED([fixed],[],[Additional fixed images for multispectral registration],[])
# ARG_OPTIONAL_REPEATED([moving],[],[Additional moving images for multispectral registration],[])
# ARG_TYPE_GROUP_SET([lineargroup],[LINEAR],[linear-type],[rigid,lsq6,similarity,lsq9,affine,lsq12,exhaustive-affine])
# ARG_OPTIONAL_SINGLE([convergence],[],[Convergence stopping value for registration],[1e-6])
# ARG_OPTIONAL_SINGLE([syn-control],[],[Non-linear (SyN) gradient and regularization parameters, not checked for correctness],[0.1,3,0])
# ARG_OPTIONAL_BOOLEAN([mask-extract],[],[Use masks to extract input images, only works with both images masked],[])
# ARG_OPTIONAL_BOOLEAN([keep-mask-after-extract],[],[Keep using masks for metric after extraction],[off])
# ARG_OPTIONAL_BOOLEAN([histogram-matching],[],[Enable histogram matching],[])
# ARG_OPTIONAL_BOOLEAN([skip-linear],[],[Skip the linear registration stages])
# ARG_OPTIONAL_BOOLEAN([skip-nonlinear],[],[Skip the nonlinear stage])
# ARG_OPTIONAL_BOOLEAN([fast],[],[Run fast SyN registration])
# ARG_OPTIONAL_BOOLEAN([float],[],[Calculate registration using float instead of double])
# ARG_OPTIONAL_BOOLEAN([clobber],[c],[Overwrite files that already exist])
# ARG_OPTIONAL_BOOLEAN([verbose],[v],[Run commands verbosely],[on])
# ARG_OPTIONAL_BOOLEAN([debug],[d],[Show all internal comands and logic for debug],[])
# ARGBASH_SET_INDENT([  ])
# ARGBASH_GO()
# needed because of Argbash --> m4_ignore([
### START OF CODE GENERATED BY Argbash v2.10.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info


die()
{
  local _ret="${2:-1}"
  test "${_PRINT_HELP:-no}" = yes && print_help >&2
  echo "$1" >&2
  exit "${_ret}"
}

# validators

lineargroup()
{
	local _allowed=("rigid" "lsq6" "similarity" "lsq9" "affine" "lsq12" "exhaustive-affine") _seeking="$1"
	for element in "${_allowed[@]}"
	do
		test "$element" = "$_seeking" && echo "$element" && return 0
	done
	die "Value '$_seeking' (of argument '$2') doesn't match the list of allowed values: 'rigid', 'lsq6', 'similarity', 'lsq9', 'affine', 'lsq12' and 'exhaustive-affine'" 4
}


begins_with_short_option()
{
  local first_option all_short_options='hocvd'
  first_option="${1:0:1}"
  test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - POSITIONALS
_positionals=()
# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_moving_mask="NOMASK"
_arg_fixed_mask="NOMASK"
_arg_resampled_output=
_arg_resampled_linear_output=
_arg_initial_transform="NONE"
_arg_linear_type="affine"
_arg_close="off"
_arg_fixed=()
_arg_moving=()
_arg_convergence="1e-6"
_arg_syn_control="0.1,3,0"
_arg_mask_extract="off"
_arg_keep_mask_after_extract="off"
_arg_histogram_matching="off"
_arg_skip_linear="off"
_arg_skip_nonlinear="off"
_arg_fast="off"
_arg_float="off"
_arg_clobber="off"
_arg_verbose="on"
_arg_debug="off"


print_help()
{
  printf '%s\n' "The general script's help msg"
  printf 'Usage: %s [-h|--help] [--moving-mask <arg>] [--fixed-mask <arg>] [-o|--resampled-output <arg>] [--resampled-linear-output <arg>] [--initial-transform <arg>] [--linear-type <LINEAR>] [--(no-)close] [--fixed <arg>] [--moving <arg>] [--convergence <arg>] [--syn-control <arg>] [--(no-)mask-extract] [--(no-)keep-mask-after-extract] [--(no-)histogram-matching] [--(no-)skip-linear] [--(no-)skip-nonlinear] [--(no-)fast] [--(no-)float] [-c|--(no-)clobber] [-v|--(no-)verbose] [-d|--(no-)debug] <movingfile> <fixedfile> <outputbasename>\n' "$0"
  printf '\t%s\n' "<movingfile>: The moving image"
  printf '\t%s\n' "<fixedfile>: The fixed image"
  printf '\t%s\n' "<outputbasename>: The basename for the output transforms"
  printf '\t%s\n' "-h, --help: Prints help"
  printf '\t%s\n' "--moving-mask: Mask for moving image (default: 'NOMASK')"
  printf '\t%s\n' "--fixed-mask: Mask for fixed image (default: 'NOMASK')"
  printf '\t%s\n' "-o, --resampled-output: Output resampled file (no default)"
  printf '\t%s\n' "--resampled-linear-output: Output resampled file with only linear transform (no default)"
  printf '\t%s\n' "--initial-transform: Initial linear transform (default: 'NONE')"
  printf '\t%s\n' "--linear-type: Type of affine transform. Can be one of: 'rigid', 'lsq6', 'similarity', 'lsq9', 'affine', 'lsq12' and 'exhaustive-affine' (default: 'affine')"
  printf '\t%s\n' "--close, --no-close: Images are starting off close, skip large scale pyramid search (off by default)"
  printf '\t%s\n' "--fixed: Additional fixed images for multispectral registration (empty by default)"
  printf '\t%s\n' "--moving: Additional moving images for multispectral registration (empty by default)"
  printf '\t%s\n' "--convergence: Convergence stopping value for registration (default: '1e-6')"
  printf '\t%s\n' "--syn-control: Non-linear (SyN) gradient and regularization parameters, not checked for correctness (default: '0.1,3,0')"
  printf '\t%s\n' "--mask-extract, --no-mask-extract: Use masks to extract input images, only works with both images masked (off by default)"
  printf '\t%s\n' "--keep-mask-after-extract, --no-keep-mask-after-extract: Keep using masks for metric after extraction (off by default)"
  printf '\t%s\n' "--histogram-matching, --no-histogram-matching: Enable histogram matching (off by default)"
  printf '\t%s\n' "--skip-linear, --no-skip-linear: Skip the linear registration stages (off by default)"
  printf '\t%s\n' "--skip-nonlinear, --no-skip-nonlinear: Skip the nonlinear stage (off by default)"
  printf '\t%s\n' "--fast, --no-fast: Run fast SyN registration (off by default)"
  printf '\t%s\n' "--float, --no-float: Calculate registration using float instead of double (off by default)"
  printf '\t%s\n' "-c, --clobber, --no-clobber: Overwrite files that already exist (off by default)"
  printf '\t%s\n' "-v, --verbose, --no-verbose: Run commands verbosely (on by default)"
  printf '\t%s\n' "-d, --debug, --no-debug: Show all internal comands and logic for debug (off by default)"
}


parse_commandline()
{
  _positionals_count=0
  while test $# -gt 0
  do
    _key="$1"
    case "$_key" in
      -h|--help)
        print_help
        exit 0
        ;;
      -h*)
        print_help
        exit 0
        ;;
      --moving-mask)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_moving_mask="$2"
        shift
        ;;
      --moving-mask=*)
        _arg_moving_mask="${_key##--moving-mask=}"
        ;;
      --fixed-mask)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_fixed_mask="$2"
        shift
        ;;
      --fixed-mask=*)
        _arg_fixed_mask="${_key##--fixed-mask=}"
        ;;
      -o|--resampled-output)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_resampled_output="$2"
        shift
        ;;
      --resampled-output=*)
        _arg_resampled_output="${_key##--resampled-output=}"
        ;;
      -o*)
        _arg_resampled_output="${_key##-o}"
        ;;
      --resampled-linear-output)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_resampled_linear_output="$2"
        shift
        ;;
      --resampled-linear-output=*)
        _arg_resampled_linear_output="${_key##--resampled-linear-output=}"
        ;;
      --initial-transform)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_initial_transform="$2"
        shift
        ;;
      --initial-transform=*)
        _arg_initial_transform="${_key##--initial-transform=}"
        ;;
      --linear-type)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_linear_type="$(lineargroup "$2" "linear-type")" || exit 1
        shift
        ;;
      --linear-type=*)
        _arg_linear_type="$(lineargroup "${_key##--linear-type=}" "linear-type")" || exit 1
        ;;
      --no-close|--close)
        _arg_close="on"
        test "${1:0:5}" = "--no-" && _arg_close="off"
        ;;
      --fixed)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_fixed+=("$2")
        shift
        ;;
      --fixed=*)
        _arg_fixed+=("${_key##--fixed=}")
        ;;
      --moving)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_moving+=("$2")
        shift
        ;;
      --moving=*)
        _arg_moving+=("${_key##--moving=}")
        ;;
      --convergence)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_convergence="$2"
        shift
        ;;
      --convergence=*)
        _arg_convergence="${_key##--convergence=}"
        ;;
      --syn-control)
        test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
        _arg_syn_control="$2"
        shift
        ;;
      --syn-control=*)
        _arg_syn_control="${_key##--syn-control=}"
        ;;
      --no-mask-extract|--mask-extract)
        _arg_mask_extract="on"
        test "${1:0:5}" = "--no-" && _arg_mask_extract="off"
        ;;
      --no-keep-mask-after-extract|--keep-mask-after-extract)
        _arg_keep_mask_after_extract="on"
        test "${1:0:5}" = "--no-" && _arg_keep_mask_after_extract="off"
        ;;
      --no-histogram-matching|--histogram-matching)
        _arg_histogram_matching="on"
        test "${1:0:5}" = "--no-" && _arg_histogram_matching="off"
        ;;
      --no-skip-linear|--skip-linear)
        _arg_skip_linear="on"
        test "${1:0:5}" = "--no-" && _arg_skip_linear="off"
        ;;
      --no-skip-nonlinear|--skip-nonlinear)
        _arg_skip_nonlinear="on"
        test "${1:0:5}" = "--no-" && _arg_skip_nonlinear="off"
        ;;
      --no-fast|--fast)
        _arg_fast="on"
        test "${1:0:5}" = "--no-" && _arg_fast="off"
        ;;
      --no-float|--float)
        _arg_float="on"
        test "${1:0:5}" = "--no-" && _arg_float="off"
        ;;
      -c|--no-clobber|--clobber)
        _arg_clobber="on"
        test "${1:0:5}" = "--no-" && _arg_clobber="off"
        ;;
      -c*)
        _arg_clobber="on"
        _next="${_key##-c}"
        if test -n "$_next" -a "$_next" != "$_key"
        then
          { begins_with_short_option "$_next" && shift && set -- "-c" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
        fi
        ;;
      -v|--no-verbose|--verbose)
        _arg_verbose="on"
        test "${1:0:5}" = "--no-" && _arg_verbose="off"
        ;;
      -v*)
        _arg_verbose="on"
        _next="${_key##-v}"
        if test -n "$_next" -a "$_next" != "$_key"
        then
          { begins_with_short_option "$_next" && shift && set -- "-v" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
        fi
        ;;
      -d|--no-debug|--debug)
        _arg_debug="on"
        test "${1:0:5}" = "--no-" && _arg_debug="off"
        ;;
      -d*)
        _arg_debug="on"
        _next="${_key##-d}"
        if test -n "$_next" -a "$_next" != "$_key"
        then
          { begins_with_short_option "$_next" && shift && set -- "-d" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
        fi
        ;;
      *)
        _last_positional="$1"
        _positionals+=("$_last_positional")
        _positionals_count=$((_positionals_count + 1))
        ;;
    esac
    shift
  done
}


handle_passed_args_count()
{
  local _required_args_string="'movingfile', 'fixedfile' and 'outputbasename'"
  test "${_positionals_count}" -ge 3 || _PRINT_HELP=yes die "FATAL ERROR: Not enough positional arguments - we require exactly 3 (namely: $_required_args_string), but got only ${_positionals_count}." 1
  test "${_positionals_count}" -le 3 || _PRINT_HELP=yes die "FATAL ERROR: There were spurious positional arguments --- we expect exactly 3 (namely: $_required_args_string), but got ${_positionals_count} (the last one was: '${_last_positional}')." 1
}


assign_positional_args()
{
  local _positional_name _shift_for=$1
  _positional_names="_arg_movingfile _arg_fixedfile _arg_outputbasename "

  shift "$_shift_for"
  for _positional_name in ${_positional_names}
  do
    test $# -gt 0 || break
    eval "$_positional_name=\${1}" || die "Error during argument parsing, possibly an Argbash bug." 1
    shift
  done
}

parse_commandline "$@"
handle_passed_args_count
assign_positional_args 1 "${_positionals[@]}"

# OTHER STUFF GENERATED BY Argbash
# Validation of values


### END OF CODE GENERATED BY Argbash (sortof) ### ])
# [ <-- needed because of Argbash

set -euo pipefail

if [[ ${_arg_debug} == "on" ]]; then
  set -x
fi

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS_PER_COMMAND:-$(nproc)}

tmpdir=$(mktemp -d)


#Setup exit trap for cleanup, don't do if debug
function finish() {
    if [[ ${_arg_debug} == "off" ]]; then
        rm -rf "${tmpdir}"
    fi

}
trap finish EXIT

#Add handler for failure to show where things went wrong
failure() {
    local lineno=$1
    local msg=$2
    echo "Failed at $lineno: $msg"
}
trap 'failure ${LINENO} "$BASH_COMMAND"' ERR


#Input checking
if [[ "${_arg_clobber}" == "off" ]]; then
  for file in ${_arg_outputbasename}0_GenericAffine.xfm ${_arg_outputbasename}0GenericAffine.mat \
              ${_arg_outputbasename}1_NL.xfm ${_arg_outputbasename}1Warp.nii.gz \
              ${_arg_resampled_output} ${_arg_resampled_linear_output}; do
    if [[ -s "${file}" ]]; then
      die "File ${file} already exists and --clobber not specified!"
    fi
  done
fi

if [[ ! ${#_arg_fixed[@]} -eq ${#_arg_moving[@]} ]]; then
  echo "Number of multispectral moving and fixed inputs not equal"
  echo "Got fixed=(${_arg_fixed[@]})"
  echo "Got moving=(${_arg_moving[@]})"
  exit 1
fi


#Check for minc or nifti, make appropriate adjustments of transforms
if [[ "${_arg_movingfile}" =~ .*"mnc" || "${_arg_fixedfile}" =~ .*"mnc" ]]; then
  minc_mode="--minc"
  second_stage_initial="${_arg_outputbasename}0_GenericAffine.xfm"
  second_stage_final="${_arg_outputbasename}1_NL.xfm"
  intermediate_resample="${tmpdir}/resample.mnc"
else
  minc_mode=""
  second_stage_initial="${_arg_outputbasename}0GenericAffine.mat"
  second_stage_final="${_arg_outputbasename}1Warp.nii.gz"
  intermediate_resample="${tmpdir}/resample.h5"
fi

#Enable verbosity
if [[ ${_arg_verbose} == "on" ]]; then
  _arg_verbose="--verbose"
else
  _arg_verbose=""
fi

#Enable histogram matching
if [[ ${_arg_histogram_matching} == "on" ]]; then
  _arg_histogram_matching=1
else
  _arg_histogram_matching=0
fi

#Float mode switch for antsRegistration
if [[ ${_arg_float} == "on" ]]; then
  _arg_float="--float 1"
else
  _arg_float="--float 0"
fi

if [[ ${_arg_mask_extract} == "on" && ${_arg_fixed_mask} != "NOMASK" && ${_arg_moving_mask} != "NOMASK" ]]; then
  ImageMath 3 ${tmpdir}/fixed_extracted.h5 m ${_arg_fixedfile} ${_arg_fixed_mask}
  ImageMath 3 ${tmpdir}/moving_extracted.h5 m ${_arg_movingfile} ${_arg_moving_mask}
  movingfile1=${tmpdir}/moving_extracted.h5
  fixedfile1=${tmpdir}/fixed_extracted.h5
  if [[ ${_arg_keep_mask_after_extract} = "off" ]]; then
    movingmask=NOMASK
    fixedmask=NOMASK
  else
    movingmask=${_arg_moving_mask}
    fixedmask=${_arg_fixed_mask}
  fi
else
  movingfile1=${_arg_movingfile}
  fixedfile1=${_arg_fixedfile}
  i=0
  while (( ${i} < ${#_arg_fixed[@]} )); do
    declare "movingfile$((i+2))=${_arg_moving[${i}]}"
    declare "fixedfile$((i+2))=${_arg_fixed[${i}]}"
    ((++i))
  done
  movingmask=${_arg_moving_mask}
  fixedmask=${_arg_fixed_mask}
fi

if [[ ${fixedmask} == "NOMASK" && ${movingmask} == "NOMASK" ]]; then
  _no_masks="--no-masks"
fi

fixed_minimum_resolution=$(python -c "print(min([abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile1} 1)\".split(\"x\")]]))")

#Calculate Maximum FOV using the size of the fixed image
#fixed_maximum_resolution=$(python -c "print(max([ a*b for a,b in zip([abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile} 1)\".split(\"x\")]],[abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile} 2)\".split(\"x\")]])]))")

#Calculate Maximum FOV using the foreground/background of the fixed image
ThresholdImage 3 ${fixedfile1} ${tmpdir}/bgmask.h5 1e-12 Inf 1 0
ThresholdImage 3 ${fixedfile1} ${tmpdir}/otsu.h5 Otsu 4 ${tmpdir}/bgmask.h5
ThresholdImage 3 ${tmpdir}/otsu.h5 ${tmpdir}/otsu.h5 2 Inf 1 0
LabelGeometryMeasures 3 ${tmpdir}/otsu.h5 none ${tmpdir}/geometry.csv
fixed_maximum_resolution=$(python -c "print(max([ a*b for a,b in zip( [ a-b for a,b in zip( [float(x) for x in \"$(tail -1 ${tmpdir}/geometry.csv | cut -d, -f 14,16,18)\".split(\",\") ],[float(x) for x in \"$(tail -1 ${tmpdir}/geometry.csv | cut -d, -f 13,15,17)\".split(\",\") ])],[abs(x) for x in [float(x) for x in \"$(PrintHeader ${fixedfile1} 1)\".split(\"x\")]])]))")

if [[ ${_arg_close} == "on" ]]; then
  steps_affine=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --convergence ${_arg_convergence} --output ${_arg_linear_type} --close ${_no_masks:+--no-masks}  --reg-pairs $((${#_arg_fixed[@]} + 1)))
else
  steps_affine=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --convergence ${_arg_convergence} --output ${_arg_linear_type} ${_no_masks:+--no-masks} --reg-pairs $((${#_arg_fixed[@]} + 1)))
fi
steps_syn=$(ants_generate_iterations.py --min ${fixed_minimum_resolution} --max ${fixed_maximum_resolution} --convergence ${_arg_convergence})

if [[ "${_arg_initial_transform}" != "NONE" ]]; then
  initial_transform=${_arg_initial_transform}
elif [[ ${_arg_skip_linear} == "off" ]]; then
  initial_transform="[ ${fixedfile1},${movingfile1},1 ]"
fi

if [[ ${_arg_skip_linear} == "off" ]]; then
  antsRegistration --dimensionality 3 ${_arg_verbose} ${minc_mode} ${_arg_float} \
    --output [ ${_arg_outputbasename} ] \
    --use-histogram-matching ${_arg_histogram_matching} \
    --initial-moving-transform ${initial_transform} \
    $(eval echo ${steps_affine})
else
  if [[ "${_arg_initial_transform}" != "NONE" ]]; then
    cp -f "${_arg_initial_transform}" "${second_stage_initial}" || true
  else
    if [[ -n ${minc_mode} ]]; then
      #Generate identity transform
      param2xfm "${second_stage_initial}"
    else
      ImageMath 3 "${second_stage_initial}" MakeAffineTransform 1
    fi
  fi
fi

if [[ ${_arg_fast} == "on" ]]; then
  syn_metric="--metric Mattes[ ${fixedfile1},${movingfile1},1,32,None ]"
  i=0
  while (( i < ${#_arg_fixed[@]} )); do
    syn_metric+=" --metric Mattes[ ${_arg_fixed[${i}]},${_arg_moving[${i}]},1,32,None ]"
    ((++i))
  done
else
  syn_metric="--metric CC[ ${fixedfile1},${movingfile1},1,4,None ]"
  i=0
  while (( i < ${#_arg_fixed[@]} )); do
    syn_metric+=" --metric CC[ ${_arg_fixed[${i}]},${_arg_moving[${i}]},1,4,None ]"
    ((++i))
  done
fi

# If requested, do linear resample
if [[ ${_arg_resampled_linear_output} && ${_arg_skip_nonlinear} == "off" ]]; then
  antsApplyTransforms -d 3 ${_arg_float} -i ${_arg_movingfile} -r ${_arg_fixedfile} -t "${second_stage_initial}" -o "${intermediate_resample}" -n BSpline[5] ${_arg_verbose}
  ThresholdImage 3 "${intermediate_resample}" "${tmpdir}/clampmask.h5" 1e-12 Inf 1 0
  ImageMath 3 "${_arg_resampled_linear_output}" m "${intermediate_resample}" "${tmpdir}/clampmask.h5"
fi

if [[ ${_arg_skip_nonlinear} == "off" ]]; then
  antsRegistration --dimensionality 3 ${_arg_verbose} ${minc_mode} ${_arg_float} \
    --output [ ${_arg_outputbasename} ] \
    --use-histogram-matching ${_arg_histogram_matching} \
    --initial-moving-transform "${second_stage_initial}" \
    --transform SyN[ ${_arg_syn_control} ] \
    ${syn_metric} \
    $(eval echo ${steps_syn}) \
    --masks [ ${fixedmask},${movingmask} ]
fi

if [[ ${_arg_resampled_output} ]]; then
  if [[ ${_arg_skip_nonlinear} == "off" ]]; then
    antsApplyTransforms -d 3 ${_arg_float} -i ${_arg_movingfile} -r ${_arg_fixedfile} -t "${second_stage_final}" -t "${second_stage_initial}" -o "${intermediate_resample}" -n BSpline[5] ${_arg_verbose}
  else
    antsApplyTransforms -d 3 ${_arg_float} -i ${_arg_movingfile} -r ${_arg_fixedfile} -t "${second_stage_initial}" -o "${intermediate_resample}" -n BSpline[5] ${_arg_verbose}
  fi
  ThresholdImage 3 "${intermediate_resample}" "${tmpdir}/clampmask.h5" 1e-12 Inf 1 0
  ImageMath 3 "${_arg_resampled_output}" m "${intermediate_resample}" "${tmpdir}/clampmask.h5"
fi

# ] <-- needed because of Argbash
