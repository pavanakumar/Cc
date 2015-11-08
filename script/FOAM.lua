-- This is the OpenFOAM reader

local ffi = require("ffi")

ffi.cdef[[
  void init_of_mesh( int *ipar );
]]

local libFOAM = ffi.load("libOpenFOAM.dylib", true)
libReader   = ffi.load("libCc_ofreader.dylib", true)

local ipar = ffi.new("int[1]", {1})
libReader.init_of_mesh( ipar )

