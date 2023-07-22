const { Polynomial } = require("../src/polynomial/polynomial.js");

function getRandomValue(higher = 10) {
    if(higher < 1) return 1;
    return Math.ceil(Math.random() * higher);
}

function getRandomArray(length, curve) {
    let buffer = [];
    for (let i = 0; i < length; i++) {
        buffer[i] = curve.Fr.random();
    }
    return buffer;
}

function getRandomBuffer(length, curve) {
    let buffer = new Uint8Array(length * curve.Fr.n8);
    for (let i = 0; i < length; i++) {
        buffer.set(curve.Fr.random(), i * curve.Fr.n8);
    }
    return buffer;
}

function getRandomPolynomialByLength(degree, curve) {
    return new Polynomial(getRandomBuffer(2 ** degree, curve), curve);
}

module.exports.getRandomValue = getRandomValue;
module.exports.getRandomArray = getRandomArray;
module.exports.getRandomBuffer = getRandomBuffer;
module.exports.getRandomPolynomialByLength = getRandomPolynomialByLength;
