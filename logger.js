const Logger = require("logplease");
const logger = Logger.create("", { showTimestamp: false });
Logger.setLogLevel("INFO");

module.exports = logger;
