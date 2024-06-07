{ pkgs ? import <nixpkgs> {}}:
let
  fhs = pkgs.buildFHSUserEnv {
    name = "crystacean";

    targetPkgs = _: [
      pkgs.micromamba
    ];

    profile = ''
      set -e
      eval "$(micromamba shell hook --shell=posix)"
      export MAMBA_ROOT_PREFIX=${builtins.getEnv "PWD"}/.mamba
      if ! test -d $MAMBA_ROOT_PREFIX/envs/mamba-crystacean; then
          micromamba create --yes -q -n mamba-crystacean
      fi
      micromamba activate mamba-crystacean
      micromamba install --yes -f requirements.txt -c conda-forge
      set +e
    '';
  };
in fhs.env