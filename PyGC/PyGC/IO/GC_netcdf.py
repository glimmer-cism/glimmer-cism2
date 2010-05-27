# Copyright 2004, Magnus Hagdorn
# 
# This file is part of GLIMMER.
# 
# GLIMMER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GLIMMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GLIMMER; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

__all__ = ['NetCDFFile']

HAVE_NCFILE=False

if not HAVE_NCFILE:
    HAVE_NCFILE=True
    try:
        from Scientific.IO.NetCDF import NetCDFFile
    except:
        HAVE_NCFILE=False

if not HAVE_NCFILE:
    HAVE_NCFILE=True
    try:
        from netCDF4 import Dataset as NetCDFFile
    except:
        HAVE_NCFILE=False

if not HAVE_NCFILE:
    raise ImportError, "Require either netCDF4 or Scientific.IO.NetCDF module"
