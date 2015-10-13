lut_file <- system.file("extdata", "color-lut.csv", package="pixelgram")
pixmap_colors <-read.csv(lut_file, header=T, stringsAsFactors=F, strip.white=T)
devtools::use_data(pixmap_colors, internal=T, overwrite=T)


### These options are from weblogo-3.4:
#  Copyright (c) 2003-2005 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks

#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.

### Popular color codings for nucleic and amino acids. 

# bp
# "TAU",  "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"
# "GC",    "blue", "Strong (3 Watson-Crick hydrogen bonds)"

# hydrophobicity
# "RKDENQ",   "blue", "hydrophilic"
# "SGHTAP",   "green", "neutral"  
# "O",   "magenta", "PNGsite"   
# "-*#",     "grey",   "unknown" 
# ".",     "white",   "transmitted" 
# "YVMCLFIW", "black",  "hydrophobic"

# chemistry
# "GSTYC",  "green",   "polar"
# "NQ",      "purple", "neutral" 
# "KRH",     "blue",   "basic"
# "O",     "cyan",   "PNGsite"
# "-*#",     "grey",   "unknown" 
# ".",     "white",   "transmitted" 
# "DE",      "red",    "acidic"
#"PAWFLIMV", "black",  "hydrophobic"

# taylor "W.R. Taylor, Protein Engineering, Vol 10, 743-746 (1997)" - check out this paper!
#   "a", "#CCFF00" 
#   "c", "#FFFF00" 
#   "d", "#FF0000"
#   "e", "#FF0066" 
#   "f", "#00FF66"
#   "g", "#FF9900"
#   "h", "#0066FF"
#   "i", "#66FF00"
#   "k", "#6600FF"
#   "l", "#33FF00"
#   "m", "#00FF00"
#   "n", "#CC00FF"
#   "o", "#00FFFF" 
#   "p", "#FFCC00"
#   "q", "#FF00CC"
#   "r", "#0000FF"
#   "s", "#FF3300"
#   "t", "#FF6600"
#   "v", "#99FF00"
#   "w", "#00CCFF"
#   "x", "#FF00FF"
#   "y", "#00FFCC"
#   "z", "#FF00FF"
#   "-", "grey" 
#   "*", "grey" 
#   "#", "grey" 
#   ".", "white" 
