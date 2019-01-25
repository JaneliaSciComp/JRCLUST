using JSON

function field2str(val)
    if isa(val, String)
        retval = "'$val'"
    elseif isa(val, Array{Any})
        retval = replace("$val", "Any" => "")
    else
        retval = "$val"
    end

    retval
end

basedir = abspath(joinpath(@__DIR__, "..")) # /*/JRCLUST

params = JSON.parse(read(joinpath(basedir, "params.json"), String))
old2new = JSON.parse(read(joinpath(basedir, "old2new.json"), String))
new2old = Dict([old2new[v] => v for v in keys(old2new)])

sections = ["usage", "probe", "recording file", "preprocessing", "spike detection", "feature extraction", "clustering", "display", "trial"]

open(joinpath(basedir, "default.prm"), "w") do io
    println(io, "% Default parameters for JRCLUST\n")
    for section in sections
        println(io, "% $(uppercase(section)) PARAMETERS")
        for (param, pdata) in params["common parameters"]
            if pdata["section"][1] == section
                print(io, "$param = $(field2str(pdata["default_value"])); % ")
                if haskey(new2old, param)
                    print(io, "(formerly $(new2old[param])) ")
                end
                if isempty(pdata["comment"])
                    println(io, pdata["description"])
                else
                    println(io, "$(pdata["description"]) ($(pdata["comment"]))")
                end
            end
        end
        if section !== "trial"
            println(io, "")
        end
    end
end
