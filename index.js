const yargs = require("yargs");

const kzg_basic_prover = require("./src/kzg_basic_prover.js");
const kzg_basic_verifier = require("./src/kzg_basic_verifier.js");

yargs
    .scriptName("kzg")
    .version("1.1.0")
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
    return await kzg_basic_prover();
}

async function _kzg_basic_verifier(proof) {
    return await kzg_basic_verifier(proof);
}
