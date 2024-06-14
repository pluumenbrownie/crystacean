{ pkgs ? import <nixpkgs> {} }:
  let
    overrides = (builtins.fromTOML (builtins.readFile ./rust-toolchain.toml));
    libPath = with pkgs; lib.makeLibraryPath [
      # load external libraries that you need in your rust project here
    ];
in
  pkgs.mkShell rec {
    buildInputs = with pkgs; [
      clang
      # Replace llvmPackages with llvmPackages_X, where X is the latest LLVM version (at the time of writing, 16)
      llvmPackages.bintools
      rustup
    ];
    RUSTC_VERSION = overrides.toolchain.channel;
    # https://github.com/rust-lang/rust-bindgen#environment-variables
    LIBCLANG_PATH = pkgs.lib.makeLibraryPath [ pkgs.llvmPackages_latest.libclang.lib ];
    shellHook = ''
      export PATH=$PATH:''${CARGO_HOME:-~/.cargo}/bin
      export PATH=$PATH:''${RUSTUP_HOME:-~/.rustup}/toolchains/$RUSTC_VERSION-x86_64-unknown-linux-gnu/bin/
      '';
    # Add precompiled library to rustc search path
    RUSTFLAGS = (builtins.map (a: ''-L ${a}/lib'') [
      # add libraries here (e.g. pkgs.libvmi)
    ]);
    LD_LIBRARY_PATH = libPath;
    # Add glibc, clang, glib, and other headers to bindgen search path
    BINDGEN_EXTRA_CLANG_ARGS =
    # Includes normal include path
    (builtins.map (a: ''-I"${a}/include"'') [
      # add dev libraries here (e.g. pkgs.libvmi.dev)
      pkgs.glibc.dev
    ])
    # Includes with special directory paths
    ++ [
      ''-I"${pkgs.llvmPackages_latest.libclang.lib}/lib/clang/${pkgs.llvmPackages_latest.libclang.version}/include"''
      ''-I"${pkgs.glib.dev}/include/glib-2.0"''
      ''-I${pkgs.glib.out}/lib/glib-2.0/include/''
    ];

    packages = [
    (pkgs.python312.withPackages (python-pkgs: [
      # select Python packages here
      # python-pkgs.maturin
      python-pkgs.typer
      python-pkgs.scipy
      python-pkgs.matplotlib
      python-pkgs.tqdm
      python-pkgs.alive-progress
      python-pkgs.ase
    ]))
    pkgs.maturin
  ];
  }
# { pkgs ? import <nixpkgs> {}}:
# let
#   fhs = pkgs.buildFHSUserEnv {
#     name = "crystacean";

#     targetPkgs = _: [
#       pkgs.micromamba
#       pkgs.cargo
#       pkgs.rustc
#       pkgs.clippy 
#       pkgs.cargo-flamegraph
#       pkgs.rustfmt
#       pkgs.gccgo
#     ];

#     profile = ''
#       export RUST_SRC_PATH="${pkgs.rust.packages.stable.rustPlatform.rustLibSrc}/lib/rustlib/src/rust/library/"
#       set -e
#       eval "$(micromamba shell hook --shell=posix)"
#       export MAMBA_ROOT_PREFIX=${builtins.getEnv "PWD"}/.mamba
#       if ! test -d $MAMBA_ROOT_PREFIX/envs/mamba-crystacean; then
#           micromamba create --yes -q -n mamba-crystacean
#       fi
#       micromamba activate mamba-crystacean
#       micromamba install --yes -f requirements.txt -c conda-forge
#       set +e
#     '';
#   };
# in fhs.env