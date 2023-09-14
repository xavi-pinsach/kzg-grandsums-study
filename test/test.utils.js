const { Polynomial } = require("../src/polynomial/polynomial.js");

function getRandomValue(min, max) {
    if(max < 1) return min;

    const random = Math.ceil(Math.random() * max);
    return Math.max(min, random);
}

function getRandomArray(length, curve) {
    let array = [];
    for (let i = 0; i < length; i++) {
        array[i] = curve.Fr.random();
    }
    return array;
}

function getBufferFromArray(array, curve) {
    let buffer = new Uint8Array(array.length * curve.Fr.n8);
    for (let i = 0; i < array.length; i++) {
        buffer.set(array[i], i * curve.Fr.n8);
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

function getRandomBinBuffer(length, curve) {
    let buffer = new Uint8Array(length * curve.Fr.n8);
    for (let i = 0; i < length; i++) {
        const bit = Math.floor(Math.random() * 2);
        buffer.set(bit === 1 ? curve.Fr.one : curve.Fr.zero, i * curve.Fr.n8);
    }
    return buffer;
}

function getOneBuffer(length, curve) {
    let buffer = new Uint8Array(length * curve.Fr.n8);
    for (let i = 0; i < length; i++) {
        buffer.set(curve.Fr.one, i * curve.Fr.n8);
    }
    return buffer;
}

function getRandomPolynomialByLength(degree, curve) {
    return new Polynomial(getRandomBuffer(2 ** degree, curve), curve);
}

module.exports.getRandomValue = getRandomValue;
module.exports.getRandomArray = getRandomArray;
module.exports.getBufferFromArray = getBufferFromArray;
module.exports.getRandomBuffer = getRandomBuffer;
module.exports.getRandomBinBuffer = getRandomBinBuffer;
module.exports.getOneBuffer = getOneBuffer;
module.exports.getRandomPolynomialByLength = getRandomPolynomialByLength;
