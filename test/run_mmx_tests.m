clear variables
close all
clc

% jump to the location of this m file
testfolder = fileparts(which(mfilename()));

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
% import matlab.unittest.plugins.StopOnFailuresPlugin

suite = TestSuite.fromFile(fullfile(testfolder,'testmmx.m'));
runner = TestRunner.withTextOutput;
% runner.addPlugin(StopOnFailuresPlugin)
result = runner.run(suite);

