#!/usr/bin/env julia
"""
Script for generating default.prm from params.json
"""

using JSON

"""
   val2str(val::Any)

Convert a Julia value to a string containing a Matlab eval-able value.
"""
function val2str(val::Any)
   if isa(val, Array{Any})
      emptystrings = isa.(val, String) .& isempty.(val)
      val[emptystrings] .= "''"
      if any(isa.(val, String)) # swap out an array for a cell array
         val = replace(replace(string(val), "Any[" => "{"), "]" => "}")
      else
         val = replace(string(val), "Any" => "")
      end
      val = replace(val, "\"''\"" => "''") # swap out "''" for ''
   elseif isa(val, String)
      val = "'$val'"
   end

   val
end

sectionorder = ["usage", "execution", "probe", "recording file",
                "preprocessing", "spike detection",
                "feature extraction", "clustering", "curation",
                "display", "trial", "validation", "preview",
                "traces", "lfp", "aux channel"]

basedir = normpath(joinpath(dirname(@__FILE__), ".."))
paramfile = joinpath(basedir, "json", "params.json")
params = JSON.parse(read(paramfile, String))
commonbysection = Dict([section => [param for param in params["common parameters"]
                        if param[2]["section"][1] == section] for section in sectionorder])
new2old = Dict([p[2] => p[1] for p in params["old2new"]])

open(joinpath(basedir, "default.prm"), "w") do io
   println(io, "% JRCLUST parameters (default parameter set)")
   println(io, "% For a description of these parameters, including default and legal values, see https://jrclust.readthedocs.io/en/latest/parameters/index.html\n")

   for section in sectionorder
      sectionparams = Dict(commonbysection[section])
      if isempty(sectionparams)
         continue
      end

      println(io, uppercase("% $section parameters"))
      for paramname in sort(String.(keys(sectionparams)))
         paramdesc = sectionparams[paramname]

         # paramname = param[1]
         # paramdesc = param[2]
         defaultval = val2str(paramdesc["default_value"])
         description = paramdesc["description"]
         comment = paramdesc["comment"]

         print(io, "$paramname = $defaultval; % ")
         if haskey(new2old, paramname)
            print(io, "(formerly $(new2old[paramname])) ")
         end

         print(io, "$description")

         if !isempty(comment)
            println(io, " ($comment)")
         else
            println(io, "")
         end
      end

      println(io, "")
   end
end
