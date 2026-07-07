#!/usr/bin/env bash

# Source this file to make the local ABangle package and CLI available.
# It prepends the repo root to PYTHONPATH so `import abangle` resolves to this
# checkout, and prepends `bin/` to PATH so the CLI is available.

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Source this file instead of executing it:" >&2
  echo "  source ./abangle_setup.sh" >&2
  exit 1
fi

abangle_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
abangle_bin_dir="${abangle_root}/bin"

path_prepend() {
  local var_name="$1"
  local new_dir="$2"
  local current_value="${!var_name}"

  case ":${current_value}:" in
    *":${new_dir}:"*) ;;
    "") printf -v "${var_name}" '%s' "${new_dir}" ;;
    *) printf -v "${var_name}" '%s:%s' "${new_dir}" "${current_value}" ;;
  esac
  export "${var_name}"
}

path_prepend PYTHONPATH "${abangle_root}"
path_prepend PATH "${abangle_bin_dir}"

unset -f path_prepend
