const { Scalar, getCurveFromQ } = require("ffjavascript");

module.exports = async function readPTauHeader(fd, sections) {
    if (!sections[1]) throw new Error(fd.fileName + ": File has no  header");
    if (sections[1].length > 1)
        throw new Error(fd.fileName + ": File has more than one header");

    fd.pos = sections[1][0].p;
    const n8 = await fd.readULE32();
    const buff = await fd.read(n8);
    const q = Scalar.fromRprLE(buff);

    const curve = await getCurveFromQ(q);

    if (curve.F1.n64 * 8 != n8) throw new Error(fd.fileName + ": Invalid size");

    const power = await fd.readULE32();
    const ceremonyPower = await fd.readULE32();

    if (fd.pos - sections[1][0].p != sections[1][0].size)
        throw new Error("Invalid PTau header size");

    return { curve, power, ceremonyPower };
};
