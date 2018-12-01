classdef tDEMO_DoyleLQGExample < matlab.unittest.TestCase
    
    properties
        testDir = pwd;
        demoFolderName = 'DoyleLQGExample';
        tempFolder;
    end
    
    properties(TestParameter)
        fileToRun = {...
            'DoyleLQGExampleFH'};
    end
    
    methods(TestMethodSetup)
        function methodSetup(testCase)
            % Create wroking folder for test
            import matlab.unittest.fixtures.TemporaryFolderFixture;
            testCase.tempFolder = testCase.applyFixture(TemporaryFolderFixture);
            
            % Copy all the files to working folder
            copyfile(fullfile(demoroot,testCase.demoFolderName,'*'),testCase.tempFolder.Folder,'f');
            cd(testCase.tempFolder.Folder);
        end
    end
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase)
            % Clean up the workspace
            bdclose('all');
            close('all');
            
            % Go back to the testDir and delete the tempFolder by Force
            cd(testCase.testDir);
            [success,msg,msgid] = rmdir(testCase.tempFolder.Folder,'s');
            if ~success
                error(msgid,msg);
            end
        end
    end
    
    methods(Test)
        function testMethod(testCase,fileToRun)
            fprintf('\n=> Running file: %s',fileToRun);
            ME = [];
            try
                evalin('base',fileToRun);
            catch ME % Catch MATLAB Exception if any
            end
            testCase.verifyEmpty(ME,'Error while running the test');
            evalin('base','clear');
        end
    end
end