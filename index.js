const yargs = require("yargs");
const version = require("./package.json").version;

const mset_eq_kzg_grandproduct_prover = require("./src/grandproduct/mset_eq_kzg_prover.js");
const mset_eq_kzg_grandproduct_verifier = require("./src/grandproduct/mset_eq_kzg_verifier.js");
const mset_eq_kzg_grandsum_prover = require("./src/grandsum/mset_eq_kzg_prover.js");
const mset_eq_kzg_grandsum_verifier = require("./src/grandsum/mset_eq_kzg_verifier.js");

yargs
    .scriptName("kzg")
    .version(version)
    .usage("$0 <cmd>")
    .command(
        "prove",
        "generates snark proofs!",
        (yargs) => {},
        function (argv) {
            return _kzg_basic_prover();
        }
    )
    .command(
        "verify [proof]",
        "verifies snark proofs",
        (yargs) => {
            yargs.positional("proof", {
                type: "object",
                default: "{}",
                describe: "snark proof",
            });
        },
        function (argv) {
            return _kzg_basic_verifier(argv.proof);
        }
    )
    .help().argv;

async function _kzg_basic_prover() {
    return await mset_eq_kzg_grandproduct_prover();
}

async function _kzg_basic_verifier(proof) {
    return await mset_eq_kzg_grandproduct_verifier(proof);
}
