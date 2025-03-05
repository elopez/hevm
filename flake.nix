{
  description = "hevm";

  inputs = {
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    foundry.url = "github:shazow/foundry.nix/47f8ae49275eeff9bf0526d45e3c1f76723bb5d3";
    solidity = {
      url = "github:ethereum/solidity/8a97fa7a1db1ec509221ead6fea6802c684ee887";
      flake = false;
    };
    ethereum-tests = {
      url = "github:ethereum/tests/v13";
      flake = false;
    };
    forge-std = {
      url = "github:foundry-rs/forge-std";
      flake = false;
    };
    solc-pkgs = {
      url = "github:hellwolf/solc.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { nixpkgs, flake-utils, solidity, forge-std, ethereum-tests, foundry, solc-pkgs, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = (import nixpkgs {
          inherit system;
          overlays = [solc-pkgs.overlay];
          config = { allowBroken = true; };
        });
        solc = (solc-pkgs.mkDefault pkgs pkgs.solc_0_8_26);
        testDeps = [
          solc
          foundry.defaultPackage.${system}
          pkgs.go-ethereum
          pkgs.z3_4_12
          pkgs.cvc5
          pkgs.git
          pkgs.bitwuzla
        ];

        secp256k1-static = stripDylib (pkgs.secp256k1.overrideAttrs (attrs: {
          configureFlags = attrs.configureFlags ++ [ "--enable-static" ];
        }));

        hspkgs = ps :
          ps.haskell.packages.ghc98.override {
            overrides = hfinal: hprev: {
              with-utf8 = ps.haskell.lib.compose.overrideCabal (drv: {
                version = "1.1.0.0";
                src = pkgs.fetchFromGitHub {
                  owner = "serokell";
                  repo = "haskell-with-utf8";
                  rev = "cf6e31475da3d9f54439650a70170819daa35f54";
                  sha256 = "sha256-hxUiZbbcA6RvrVgGk4Vbt/rZT6wnBF3bfYbbQflzQ24=";
                };
              }) hprev.with-utf8;
            };
          };
        hlib = pkgs.haskell.lib;

        # base hevm derivation.
        # parameterized on the pkgs definition to allow use of `pkgsStatic` or `pkgs` as needed.
        hevmBase = ps :
          ps.lib.pipe
            (((hspkgs ps).callCabal2nix "hevm" ./. {
              secp256k1 = ps.secp256k1;
            }).overrideAttrs(final: prev: {
              HEVM_SOLIDITY_REPO = solidity;
              HEVM_ETHEREUM_TESTS_REPO = ethereum-tests;
              HEVM_FORGE_STD_REPO = forge-std;
              DAPP_SOLC = "${solc}/bin/solc";
            }))
            [
              (hlib.compose.overrideCabal (old: { testTarget = "test"; }))
              (hlib.compose.addTestToolDepends testDeps)
              (hlib.compose.appendConfigureFlags [ "-fci" "-O2" ])
              hlib.compose.dontHaddock
              hlib.compose.doCheck
            ];

        # "static" binary for distribution
        # on linux this is actually a real fully static binary
        # on macos this has everything except libcxx, libsystem and libiconv
        # statically linked. we can be confident that these three will always
        # be provided in a well known location by macos itself.
        hevmRedistributable = let
          grep = "${pkgs.gnugrep}/bin/grep";
          otool = "${pkgs.darwin.binutils.bintools}/bin/otool";
          install_name_tool = "${pkgs.darwin.binutils.bintools}/bin/install_name_tool";
          codesign_allocate = "${pkgs.darwin.binutils.bintools}/bin/codesign_allocate";
          codesign = "${pkgs.darwin.sigtool}/bin/codesign";
        in if pkgs.stdenv.isLinux
        then hlib.dontCheck (hevmBase pkgs.pkgsMusl.pkgsStatic)
        else pkgs.runCommand "stripNixRefs" {} ''
          mkdir -p $out/bin
          cp ${hlib.dontCheck (forceStaticDepsMacos (hevmBase pkgs))}/bin/hevm $out/bin/

          # get the list of dynamic libs from otool and tidy the output
          libs=$(${otool} -L $out/bin/hevm | tail -n +2 | sed 's/^[[:space:]]*//' | cut -d' ' -f1)

          # get the paths for libcxx and libiconv
          cxx=$(echo "$libs" | ${grep} '^/nix/store/.*/libc++\.')
          cxxabi=$(echo "$libs" | ${grep} '^/nix/store/.*/libc++abi\.')
          iconv=$(echo "$libs" | ${grep} '^/nix/store/.*/libiconv\.')

          # rewrite /nix/... library paths to point to /usr/lib
          chmod 777 $out/bin/hevm
          ${install_name_tool} -change "$cxx" /usr/lib/libc++.1.dylib $out/bin/hevm
          ${install_name_tool} -change "$cxxabi" /usr/lib/libc++abi.dylib $out/bin/hevm
          ${install_name_tool} -change "$iconv" /usr/lib/libiconv.dylib $out/bin/hevm
          # check that no nix deps remain
          nixdeps=$(${otool} -L $out/bin/hevm | tail -n +2 | { ${grep} /nix/store -c || test $? = 1; })
          if [ ! "$nixdeps" = "0" ]; then
            echo "Nix deps remain in redistributable binary!"
            exit 255
          fi
          # re-sign binary
          CODESIGN_ALLOCATE=${codesign_allocate} ${codesign} -f -s - $out/bin/hevm
          chmod 555 $out/bin/hevm
        '';


        # wrapped binary for use on systems with nix available.
        # does not statically link.
        # ensures all required runtime deps are available and on path.
        hevmWrapped =
          pkgs.symlinkJoin {
            name = "hevm";
            paths = [ (hlib.dontCheck (hevmBase pkgs)) ];
            buildInputs = [ pkgs.makeWrapper ];
            postBuild = ''
              wrapProgram $out/bin/hevm \
                --prefix PATH : "${pkgs.lib.makeBinPath [
                  pkgs.bash
                  pkgs.coreutils
                  pkgs.git
                  pkgs.z3
                  pkgs.cvc5
                  pkgs.bitwuzla
                  solc
                ]}"
            '';
          };

        # if we pass a library folder to ghc via --extra-lib-dirs that contains
        # only .a files, then ghc will link that library statically instead of
        # dynamically (even if --enable-executable-static is not passed to cabal).
        # we use this trick to force static linking of some libraries on macos.
        stripDylib = drv : pkgs.runCommand "${drv.name}-strip-dylibs" {} ''
          mkdir -p $out
          mkdir -p $out/lib
          cp -r ${drv}/* $out/
          rm -rf $out/**/*.dylib
        '';

        # ensures that all required deps will be linked statically on macos builds
        forceStaticDepsMacos = p :
          hlib.appendConfigureFlags p
            [
              "--extra-lib-dirs=${stripDylib (pkgs.gmp.override { withStatic = true; })}/lib"
              "--extra-lib-dirs=${stripDylib secp256k1-static}/lib"
              "--extra-lib-dirs=${stripDylib (pkgs.libff.override { enableStatic = true; })}/lib"
              "--extra-lib-dirs=${pkgs.zlib.static}/lib"
              "--extra-lib-dirs=${stripDylib (pkgs.libffi.overrideAttrs (_: { dontDisableStatic = true; }))}/lib"
              "--extra-lib-dirs=${stripDylib (pkgs.ncurses.override { enableStatic = true; })}/lib"
            ];


      in rec {

        # --- packages ----

        packages.ci = pkgs.lib.pipe (hevmBase (if pkgs.stdenv.isLinux then pkgs.pkgsMusl.pkgsStatic else pkgs)) (with hlib.compose; [doBenchmark dontHaddock disableLibraryProfiling]);
        packages.unwrapped = hlib.dontCheck (hevmBase pkgs);
        packages.hevm = hevmWrapped;
        packages.redistributable = hevmRedistributable;
        packages.default = packages.hevm;

        # --- apps ----

        apps.hevm = flake-utils.lib.mkApp { drv = packages.hevm; };
        apps.default = apps.hevm;

        # --- shell ---

        devShells.default = let
          libraryPath = "${pkgs.lib.makeLibraryPath [ pkgs.libff pkgs.secp256k1 pkgs.gmp ]}";
        in (hspkgs pkgs).shellFor {
          packages = _: [ (hevmBase pkgs) ];
          buildInputs = [
            pkgs.curl
            pkgs.mdbook
            pkgs.yarn
            (hspkgs pkgs).cabal-install
            (hspkgs pkgs).eventlog2html
            (hspkgs pkgs).haskell-language-server
          ] ++ testDeps;
          withHoogle = true;

          # hevm tests expect these to be set
          HEVM_SOLIDITY_REPO = solidity;
          DAPP_SOLC = "${solc}/bin/solc";
          HEVM_ETHEREUM_TESTS_REPO = ethereum-tests;
          HEVM_FORGE_STD_REPO = forge-std;

          # point cabal repl to system deps
          LD_LIBRARY_PATH = libraryPath;
          shellHook = pkgs.lib.optionalString pkgs.stdenv.isDarwin ''
            export DYLD_LIBRARY_PATH="${libraryPath}";
          '';
        };
      }
    );
}
