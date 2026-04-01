function testlog(loglevel, logfilepath, targetfile, verbose)
    if loglevel == Logging.Warn
        @testset verbose = verbose ≥ 2 "Test log" begin
            testlog = open(logfilepath, "r") do io
                read(io, String)
            end
            targetlog = open(targetfile, "r") do io
                read(io, String)
            end
            @test testlog == targetlog
        end
    end
end
