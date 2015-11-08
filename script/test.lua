--- LuaJIT program to run the reader and FORTRAN Cc code

require("FOAM")

local ffi = require("ffi")

local pm_ptr = ffi.new("void *")


