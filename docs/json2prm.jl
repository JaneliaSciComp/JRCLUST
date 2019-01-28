using JSON

function field2str(val)
    if isa(val, String)
        retval = "'$val'"
    # elseif isa(val, Array{Any}) && all([isa(v, String) for v in val])
    #     retval = "{'$(join(val, "', '"))'}"
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

sections = ["usage", "execution", "probe", "recording file", "preprocessing", "spike detection", "feature extraction", "clustering", "curation", "display", "trial", "validation", "preview", "traces", "lfp", "aux channel"]

open(joinpath(basedir, "docs", "source", "parameters", "index.rst"), "w") do io
    println(io, ".. _parameters:\n")
    println(io, "JRCLUST parameters\n==================\n")
    println(io, "Common parameters\n------------------\n")
    common = params["common parameters"]
    # for section in sections
    #     header = "$(titlecase(section)) parameters"
    #     println(io, header)
    #     println(io, "$(repeat("~", length(header)))\n")
        commonkeys = sort([k for k in keys(common)])
        for param in commonkeys
            pdata = common[param]
            # if pdata["section"][1] == section
                println(io, ".. _$param:\n")
                ptitle = "``$param``"
                println(io, ptitle)
                println(io, "$(repeat("^", length(ptitle)))\n")
                if haskey(new2old, param)
                    println(io, "(Formerly ``$(new2old[param])``)\n")
                end
                # print("$param = $(field2str(pdata["default_value"])); % ")
                if isempty(pdata["comment"])
                    println(io, "$(pdata["description"]).\n")
                else
                    println(io, "$(pdata["description"]) ($(pdata["comment"])).\n")
                end

                pvalid = pdata["validation"]
                if haskey(pvalid, "values")
                    println(io, "One of the following:\n")
                    for val in pvalid["values"]
                        println(io, "- '$val'")
                    end
                    println(io, "")
                end

                defval = field2str(pdata["default_value"])
                if defval == "''"
                    defval = "an empty string"
                elseif defval == "[]"
                    defval = "empty"
                end
                println(io, "**Default** is $(defval).\n")
            # end
        end
    # end

    println(io, "Advanced parameters\n-------------------\n")
    advanced = params["advanced parameters"]
    # for section in sections
    #     header = "$(titlecase(section)) parameters"
    #     println(io, header)
    #     println(io, "$(repeat("~", length(header)))\n")
        advancedkeys = sort([k for k in keys(advanced)])
        for param in advancedkeys
            pdata = advanced[param]
            # if pdata["section"][1] == section
                println(io, ".. _$param:\n")
                ptitle = "``$param``"
                println(io, ptitle)
                println(io, "$(repeat("^", length(ptitle)))\n")
                if haskey(new2old, param)
                    println(io, "(Formerly ``$(new2old[param])``)\n")
                end
                # print("$param = $(field2str(pdata["default_value"])); % ")
                if isempty(pdata["comment"])
                    println(io, "$(pdata["description"]).\n")
                else
                    println(io, "$(pdata["description"]) ($(pdata["comment"])).\n")
                end

                pvalid = pdata["validation"]
                if haskey(pvalid, "values")
                    println(io, "One of the following:\n")
                    for val in pvalid["values"]
                        println(io, "- '$val'")
                    end
                    println(io, "")
                end

                defval = field2str(pdata["default_value"])
                if defval == "''"
                    defval = "an empty string"
                elseif defval == "[]"
                    defval = "empty"
                end
                println(io, "**Default** is $(defval).\n")
            # end
        end
        # for param in paramkeys
        #     pdata = params[param]

        # end
        # if section !== "trial"
        #     println(io, "")
        # end
    # end
end

# do old2new getters/setters
# newvals = sort([k for k in keys(new2old)])
# open("foo.txt", "w") do io
#     for newval in newvals
#         oldval = new2old[newval]
#         if !haskey(params, newval)
#             println("not in params: $newval")
#             continue
#         end
#         pdata = params[newval]["validation"]
#         if isempty(pdata["attributes"])
#             println("no attributes: $newval")
#             continue
#         end
#         println(io, "        % $(newval)/$(oldval)")
#         println(io, "        function set.$(newval)(obj, val)")
#         println(io, "            validateattributes(val, {'$(join(pdata["classes"], "', '"))'}, {'$(join(pdata["attributes"], "', '"))'});")
#         if haskey(pdata, "postapply")
#             println(io, "            hFun = $(pdata["postapply"]);")
#             println(io, "            val = hFun(val);")
#         end
#         if haskey(pdata, "postassert")
#             println(io, "            hFun = $(pdata["postassert"]);")
#             println(io, "            assert(hFun(val));")
#         end
#         if haskey(pdata, "values")
#             println(io, "            legalVals = {'$(join(pdata["values"], "', '"))'};")
#             println(io, "            assert(ismember(val, legalVals), 'legal values are %s', strjoin(legalVals, ', '))")
#         end
#         println(io, "            obj.$(newval) = val;")
#         println(io, "        end")
#         println(io, "        function val = get.$(oldval)(obj)")
#         println(io, "            obj.logOldP('$(oldval)');")
#         println(io, "            val = obj.$(newval);")
#         println(io, "        end")
#         println(io, "        function set.$(oldval)(obj, val)")
#         println(io, "            obj.logOldP('$(oldval)');")
#         println(io, "            obj.$(newval) = val;")
#         println(io, "        end")
#         println(io, "        ")
#     end
# end
