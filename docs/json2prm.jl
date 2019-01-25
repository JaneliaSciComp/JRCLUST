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

params = merge(params["common parameters"], params["advanced parameters"])
paramkeys = sort([k for k in keys(params)])
sections = ["usage", "execution", "probe", "recording file", "preprocessing", "spike detection", "feature extraction", "clustering", "curation", "display", "trial", "validation", "preview", "traces", "lfp", "aux channel"]

#open(joinpath(basedir, "default.prm"), "w") do io
# println("% Default parameters for JRCLUST\n")
# for section in sections
#     println("% $(uppercase(section)) PARAMETERS")
#     for param in paramkeys
#         pdata = params[param]
#         if pdata["section"][1] == section
#             print("$param = $(field2str(pdata["default_value"])); % ")
#             println(pdata["description"])
#         end
#     end
#     if section !== "trial"
#         println("")
#     end
# end
#end

# do old2new getters/setters
newvals = sort([k for k in keys(new2old)])
open("foo.txt", "w") do io
    for newval in newvals
        oldval = new2old[newval]
        if !haskey(params, newval)
            println("not in params: $newval")
            continue
        end
        pdata = params[newval]["validation"]
        if isempty(pdata["attributes"])
            println("no attributes: $newval")
            continue
        end
        println(io, "        % $(newval)/$(oldval)")
        println(io, "        function set.$(newval)(obj, val)")
        println(io, "            validateattributes(val, {'$(join(pdata["classes"], "', '"))'}, {'$(join(pdata["attributes"], "', '"))'});")
        if haskey(pdata, "postapply")
            println(io, "            hFun = $(pdata["postapply"]);")
            println(io, "            val = hFun(val);")
        end
        if haskey(pdata, "postassert")
            println(io, "            hFun = $(pdata["postassert"]);")
            println(io, "            assert(hFun(val));")
        end
        if haskey(pdata, "values")
            println(io, "            legalVals = {'$(join(pdata["values"], "', '"))'};")
            println(io, "            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))")
        end
        println(io, "            obj.$(newval) = val;")
        println(io, "        end")
        println(io, "        function val = get.$(oldval)(obj)")
        println(io, "            obj.logOldP('$(oldval)');")
        println(io, "            val = obj.$(newval);")
        println(io, "        end")
        println(io, "        function set.$(oldval)(obj, val)")
        println(io, "            obj.logOldP('$(oldval)');")
        println(io, "            obj.$(newval) = val;")
        println(io, "        end")
        println(io, "        ")
    end
end
