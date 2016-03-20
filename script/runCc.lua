-- Lua code to run the Cc solver

local ffi = require("ffi")

ffi.cdef[[
  void init_cc(int *parflag, int *nlevel);
  void finalise_cc();
]]

local libCc = ffi.load("libCc.so", true)

x = 3
local parflag = ffi.new("int[1]", {1})
local nlevel  = ffi.new("int[1]", {x})

nlevel[0] = nlevel[0] + 1
-- libCc.init_cc( parflag, nlevel )
libCc.init_cc( parflag, nlevel )
libCc.finalise_cc();


