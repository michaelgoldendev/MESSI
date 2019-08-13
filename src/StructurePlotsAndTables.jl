using PyPlot
using StatsBase
using Distributions
using FastaIO
using Printf

function summarizemcmc(mcmclogfile)
  fin = open(mcmclogfile,"r")
  params = Dict{AbstractString,Array{Float64,1}}()
  header = true
  keys = AbstractString[]
  for line in readlines(fin)
    spl = split(strip(line))
    if header
      keys = spl
      for key in keys
        params[key] = Float64[]
      end
      header = false
    else
      for i=1:length(spl)
        p = get(params, keys[i], Float64[])
        push!(p, parse(Float64,spl[i]))
        params[keys[i]] = p
      end
    end
  end
  close(fin)
  for key in keys
    len = length(params[key])
    if len > 1
      data = params[key][div(len,2):end]
      #println(key,"\t",length(data),"\t", mean(data), "\t", std(data))
    end
  end
  len = length(params[keys[1]])
  #=
  if len > 2
    println("lambdaGC > lambdaAT: ", mean([lambdaGC > lambdaAT ? 1.0 : 0.0 for (lambdaAT,lambdaGC) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
    println("lambdaAT > lambdaGT: ", mean([lambdaAT > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaAT) in zip(params["lambdaGT"][div(len,2):end],params["lambdaAT"][div(len,2):end])]))
    println("lambdaGC > lambdaGT: ", mean([lambdaGC > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaGC) in zip(params["lambdaGT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
  end=#
  return params, len
end

function plotstructurebenchmarks()
  PyPlot.matplotlib[:rc]("text", usetex=true)
  PyPlot.matplotlib[:rcParams]["text.latex.unicode"] = true
  PyPlot.matplotlib[:rc]("font", family="serif", size=21)
  fig = figure("pyplot_histogram",figsize=(14,6))
  ax = fig.gca()

  #margins(x=0,y=0,tight=false)
  metrics = ["Precision", "Recall", "F1 score", "MCC", "Weighted mountain\nsimilarity", "Mountain\nsimilarity (p=1)"]
  nmetrics = length(metrics)
  colors = ["green", "#0f87bf","red"]
  methods = ["Our method", "RNAalifold", "PPFold"]
  nmethods = length(methods)















  ourmethod = Array{Float64,1}[[0.871429,0.855556,0.791667,0.87234,0.424242,0.469388,0.885714,0.901639,0.810811,0.54386,0.707317,0.830189,0.653333,0.640625,0.25,0.962264,0.663717,0.616667,0.857143,0.650052,0.677215,0.8],[0.67033,0.974684,0.76,0.854167,0.875,0.92,0.984127,0.632184,1.0,0.861111,0.852941,0.977778,0.720588,0.863158,1.0,0.85,0.986842,0.986667,0.934426,0.932936,0.896648,0.9],[0.757764,0.911243,0.77551,0.863158,0.571429,0.621622,0.932331,0.743243,0.895522,0.666667,0.773333,0.897959,0.685315,0.735426,0.4,0.902655,0.793651,0.758974,0.894118,0.766218,0.771635,0.847059],[0.799026,0.913142,0.775414,0.863125,0.64981,0.869229,0.933538,0.754774,0.900278,0.745638,0.776412,0.900844,0.686107,0.74356,0.512833,0.904267,0.820189,0.846034,0.894868,0.78281,0.782542,0.848217],[0.856172,0.958199,0.907464,0.966109,0.898204,0.826087,0.955414,0.853782,0.926316,0.902541,0.932933,0.927536,0.983661,0.862924,0.859813,0.902439,0.877076,0.909457,0.950019,0.857844,0.901629,0.910891],[0.935862,0.960994,0.990851,0.96281,0.868442,0.741615,0.982334,0.964381,0.957135,0.900679,0.973079,0.981573,0.990007,0.957374,0.976289,0.9384,0.908665,0.912356,0.992501,0.970572,0.978323,0.975636]]

  ppfoldmethod = Array{Float64,1}[[0.833333,0.888889,0.678571,0.930233,0.518519,0.533333,0.784314,0.903226,0.878788,0.507246,0.714286,0.913043,0.641791,0.637168,0.384615,0.959184,0.633929,0.591304,0.843972,0.676667,0.707263,0.9],[0.274725,0.607595,0.76,0.833333,0.875,0.96,0.634921,0.643678,0.966667,0.972222,0.735294,0.933333,0.632353,0.757895,1.0,0.783333,0.934211,0.906667,0.97541,0.907601,0.884078,0.9],[0.413223,0.721805,0.716981,0.879121,0.651163,0.685714,0.701754,0.751678,0.920635,0.666667,0.724638,0.923077,0.637037,0.692308,0.555556,0.862385,0.755319,0.715789,0.904943,0.775302,0.785847,0.9],[0.478003,0.73482,0.717784,0.880385,0.729712,0.941319,0.705344,0.762281,0.921548,0.833265,0.724353,0.923042,0.637021,0.694853,0.620066,0.866645,0.779975,0.797212,0.90724,0.787603,0.79429,0.899805],[0.611465,0.907114,0.897759,0.961801,0.934132,0.855072,0.889666,0.862097,0.947368,0.85,0.96359,0.949275,0.978004,0.875043,0.925234,0.882927,0.877229,0.918483,0.933099,0.881647,0.914646,0.948036],[0.686769,0.978488,0.979929,0.946776,0.892643,0.773913,0.976605,0.968688,0.987954,0.800823,0.991379,0.995859,0.991211,0.964409,0.986933,0.907906,0.918051,0.91182,0.980677,0.972837,0.980281,0.986022]]

  rnaalifoldmethod = Array{Float64,1}[[1.0,0.942857,1.0,0.967742,0.518519,0.547619,1.0,0.921053,1.0,0.75,0.84,1.0,0.705882,0.698246,1.0,0.987342,0.766667,0.666667,0.877049,0.76834,0.802632,0.897436],[0.318681,0.417722,0.48,0.625,0.875,0.92,0.603175,0.402299,0.7,0.583333,0.617647,0.888889,0.529412,0.698246,1.0,0.65,0.907895,0.72,0.877049,0.889717,0.851955,0.875],[0.483333,0.578947,0.648649,0.759494,0.651163,0.686567,0.752475,0.56,0.823529,0.65625,0.711864,0.941176,0.605042,0.698246,1.0,0.78392,0.831325,0.692308,0.877049,0.824586,0.826558,0.886076],[0.564164,0.62749,0.692609,0.777617,0.729712,0.938932,0.776448,0.608477,0.836454,0.763704,0.719994,0.942748,0.611281,0.69819,1.0,0.800898,0.834224,0.779392,0.876956,0.828927,0.831304,0.885927],[0.605096,0.856724,0.873786,0.906968,0.934076,0.876812,0.840764,0.744477,0.905263,0.945371,0.909848,0.963768,0.970497,0.888796,1.0,0.8,0.947739,0.960469,0.938609,0.927367,0.938935,0.930693],[0.692091,0.967024,0.952203,0.914241,0.903693,0.789493,0.9321,0.932748,0.979924,0.959729,0.970811,0.994462,0.988702,0.970155,1.0,0.861341,0.945473,0.935037,0.988362,0.982108,0.985524,0.981751]]

  nstructures = length(ourmethod[1])


  ciwidth = 10.0

  offsetx = 0.0
  xticklabelpos = Float64[]
  for metric=1:nmetrics
    x = [1.5*i for i=1:nmethods]
    y = [mean(ourmethod[metric]), mean(ppfoldmethod[metric]), mean(rnaalifoldmethod[metric])]
    #updownerror = [100*mean(ourmethod[metric])]
    errors = [yi*rand() for yi in y]
    #=
    ourmethod5 = y[1] - percentile(ourmethod[metric],ciwidth)
    ourmethod95 = percentile(ourmethod[metric],100.0 - ciwidth)  - y[1]
    ppfoldmethod5 = y[2] - percentile(ppfoldmethod[metric],ciwidth)
    ppfoldmethod95 = percentile(ppfoldmethod[metric],100.0 - ciwidth)  - y[2]
    rnaalifold5 =  y[3] - percentile(rnaalifoldmethod[metric],ciwidth)
    rnaalifold95 = percentile(rnaalifoldmethod[metric],100.0 - ciwidth) - y[3]
    =#
    ourmethod5 = std(ourmethod[metric])
    ourmethod95 = std(ourmethod[metric])
    ppfoldmethod5 = std(ppfoldmethod[metric])
    ppfoldmethod95 = std(ppfoldmethod[metric])
    rnaalifold5 = std(rnaalifoldmethod[metric])
    rnaalifold95 = std(rnaalifoldmethod[metric])
    println(ourmethod5, "\t", ourmethod95)
    println(ppfoldmethod5, "\t", ppfoldmethod95)
    println(rnaalifold5, "\t", rnaalifold95)


    for (xi, yi) in zip(x,y)
      t = string(L"$",@sprintf("%.2f", yi),L"$")
      text(xi+offsetx+0.05, 0.32, string(@sprintf("%.2f", yi)), fontsize=14, horizontalalignment="center")

      #annotate(t, xy=[xi+offsetx;0.1],textcoords="offset points", xytext=[xi+offsetx;0.1])
      #=textcoords="offset points",
      fontsize=16.0,
      ha="center",
      va="bottom")=#
    end

    updownerror = [[ourmethod5, ppfoldmethod5, rnaalifold5], [ourmethod95, ppfoldmethod95, rnaalifold95]]
    b = bar(x.+offsetx,y,1.35,color=colors,align="center",alpha=0.4, yerr=updownerror)
    push!(xticklabelpos, 3.0 + offsetx)
    offsetx += nmethods*2.0




  end

  #println(keys(leg))
  #=
  for (method, color) in zip(methods,colors)
    legend(methods, color=color)
  end=#

  #leg = legend(methods, facecolors=[[1.0,0.0,0.0,0.5],[0.0,1.0,0.0,0.5],[0.0,0.0,1.0,0.5]])


  xticks(xticklabelpos, metrics, rotation="0")
  ylabel("Mean score and standard deviation", fontsize=21)


  leg = legend(methods, loc=9, ncol=3)
  #println(keys(leg))
  patches = leg[:get_patches]()
  for (patch, col)  in zip(patches,colors)
    patch[:set_facecolor](col)
  end

  ax[:set_ylim]([0.3,1.15])

  title("Summary of secondary structure prediction benchmarks (N=$nstructures)")

  savefig("benchmarksprediction.svg")

end



function plotratios(outputprefix)
  params,len = summarizemcmc(string(outputprefix,".B1.0.M1.0.mcmc.log"))
  GCdivAT = [lambdaGC / lambdaAT for (lambdaAT,lambdaGC) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]
  ATdivGT = [lambdaAT / lambdaGT for (lambdaAT,lambdaGT) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGT"][div(len,2):end])]
  GCdivGT = [lambdaGC / lambdaGT for (lambdaGT,lambdaGC) in zip(params["lambdaGT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]

  PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
  PyPlot.matplotlib[:rcParams]["text.latex.unicode"] = true
  PyPlot.matplotlib[:rc]("font", family="serif", size=21)
  fig = nothing
  basey = 0.05
  leglabels = []
  fracs = [L"${\lambda_\mathrm{GC}}/{\lambda_\mathrm{AU}}$", L"${\lambda_\mathrm{AU}}/{\lambda_\mathrm{GU}}$", L"${\lambda_\mathrm{GC}}/{\lambda_\mathrm{GU}}$"]
  index = 1
  colors=["darkmagenta", "orange", "darkcyan"]
  for (x,t,f, xoffset) in zip((GCdivAT,ATdivGT,GCdivGT),(L"$\frac{\lambda_\mathrm{GC}}{\lambda_\mathrm{AU}}$", L"$\frac{\lambda_\mathrm{AU}}{\lambda_\mathrm{GU}}$", L"$\frac{\lambda_\mathrm{GC}}{\lambda_\mathrm{GU}}$"), ("GC_AU","AU_GU", "GC_GU"), (0.0,0.0,0.005))
    nbins = 15 # Number of bins

    fig = figure("pyplot_histogram",figsize=(8,7)) # Not strictly required
    ax = axes() # Not strictly required

    h = plt[:hist](x,nbins, normed=1,color=colors[index], alpha=0.7) # Histogram
    #ax[:set_xlim]([0.0625,16.0])
    ax[:set_xlim]([0.25,16.0])
    ax[:set_xscale]("log")
    grid("on")
    xlabel("Ratio (log scale)")
    ylabel("Posterior density")


    #v = (log10(median(x)) + 1.0) / 2.0
    md = string(L"$\textit{Md}=",@sprintf("%.1f", median(x)),L"$")
    #v = 0.5 + log10(median(x))/log10(16)/2.0

    #v = 0.5 + log10(median(x))/log10(16)/2.0
    v = 2.0/(8.0-2.0) + (1.0 - (2.0/(8.0-2.0)))*log10(median(x))/log10(16)
    println(v)
    annotate(t,
  	xy=[v+xoffset;basey],
  	xycoords="axes fraction",
    xytext=[0.0,0],
  	textcoords="offset points",
  	fontsize=30.0,
  	ha="center",
  	va="bottom")

    title("Posterior ratios of coevolution rates")

    #xticks([0.0625,0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0], ["0.0625", "0.125", "0.25", "0.5", "1.0", "2.0", "4.0", "8.0", "16.0"])
    xticks([0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0], ["0.25", "0.5", "1.0", "2.0", "4.0", "8.0", "16.0"])
    #xticks([0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0], [""0.25", "0.5", "1.0", "2.0", "4.0", "8.0", "16.0"])
    #methods = [L"$\lambda_\mathrm{GC}$", L"$\lambda_\mathrm{AU}$", L"$\lambda_\mathrm{GU}$"]
    push!(leglabels, string(fracs[index], " (",md,")"))
    legend(leglabels, loc=8, ncol=2,fontsize=19,bbox_to_anchor=(0.485, -0.475))
    #legend(leglabels, loc=8, ncol=2,fontsize=22,bbox_to_anchor=(0.5, -0.535))
    #margins(bottom="20cm")
    #l[:get_lines][1][:]
    fig[:subplots_adjust](bottom=0.30)

    #savefig(string(outputprefix,".figure.",f,".svg"))
    #close(fig)
    index += 1
  end
  axvline(x=1.0,color="black",linestyle="dotted")
  savefig(string(outputprefix,"_figure_ratios.svg"))
  close(fig)
end



function plotdistributions(outputprefix)
  params,len = summarizemcmc(string(outputprefix,".B1.0.M1.0.mcmc.log"))
  GC = params["lambdaGC"][div(len,2):end]
  AU = params["lambdaAT"][div(len,2):end]
  GU = params["lambdaGT"][div(len,2):end]

  PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
  PyPlot.matplotlib[:rcParams]["text.latex.unicode"] = true
  PyPlot.matplotlib[:rc]("font", family="serif", size=21)
  fig = nothing
  #L"$\tilde{\lambda}_\mathrm{GC}$", L"$\tilde{\lambda}_\mathrm{AU}$", L"$\tilde{\lambda}_\mathrm{GU}$"

  leglabels = []

  index = 1
  colors=["blue","crimson",  "green"]
  for (x,t,f, xoffset) in zip((GC,AU,GU),(L"$\lambda_\mathrm{GC}$", L"$\lambda_\mathrm{AU}$", L"$\lambda_\mathrm{GU}$"), ("GC_AU","AU_GU", "GC_GU"), (0.0,0.0,0.005))
    nbins = 15 # Number of bins

    fig = figure("pyplot_histogram",figsize=(8,7)) # Not strictly required
    ax = axes() # Not strictly required
    h = plt[:hist](x,nbins, normed=1,color=colors[index],alpha=0.6) # Histogram
    grid("on")
    ylabel("Posterior density")

    logscale = false
    #v = + log10(median(x))/log10(16)/2.0
    xmax = 16.0
    #xmax = 32.0
    v = (median(x) - 1.0) / (xmax-1.0)
    if logscale
      v = (log10(median(x))) / log10(xmax)
    end
    v = max(v,0.05)
    medianstr = @sprintf("%.1f", median(x))
    t = string(t)
    println(v)
    #=
    text(v+xoffset, 0.05, t,
    fontsize=22.0,
    ha="center",
    va="bottom")=#

    annotate(t,
  	xy=[v+xoffset;0.09],
  	xycoords="axes fraction",
    xytext=[0.0,0],
  	textcoords="offset points",
  	fontsize=26.0,
  	ha="center",
  	va="bottom")


    md = string(L"$\textit{Md}=",medianstr,L"$")
    #=
    annotate(md,
  	xy=[v+xoffset;0.04],
  	xycoords="axes fraction",
    xytext=[0.0,0],
  	textcoords="offset points",
  	fontsize=16.0,
  	ha="center",
  	va="bottom")=#

    title("Posterior distributions of coevolution rates")


    if logscale
      xlabel("Rate (log scale)")
      ax[:set_xlim]([1.0,16.0])
      ax[:set_xscale]("log")
      xticks([1.0, 2.0, 4.0, 8.0], ["1.0", "2.0", "4.0", "8.0", "16.0"])
    else
      xlabel("Rate")
      ax[:set_xlim]([1.0,xmax])
      xtickvals = [1.0]
      while xtickvals[end] < xmax
        push!(xtickvals, xtickvals[end]+3.0)
      end
      xtickstrs = [string(v) for v in xtickvals]

      xticks(xtickvals, xtickstrs)
    end

    push!(leglabels, string(t, " (",md,")"))
    legend(leglabels, loc=8, ncol=2,fontsize=19,bbox_to_anchor=(0.495, -0.475))
    #margins(bottom="20cm")
    fig[:subplots_adjust](bottom=0.30)
    index += 1

    #ax[:set_ylim]([0.0,2.0])
    println("MY YLIM=", string(ax[:get_ylim]))
  end
  savefig(string(outputprefix,"_figure_distributions.svg"))
  close(fig)
end

function printtables(files)
  ret = ""
  for maxfile in files
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]

    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    #maxZentropy = readentropy(string(maxfile, ".entropy"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZguparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))
    delta = "-"
    ret = string(ret, " & Unpaired & ", @sprintf("%.2f", maxZunpaired), " & ",delta," & 10 & ", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    delta = @sprintf("%.2f", maxZgu-maxZunpaired)
    pval = 1.0 - cdf(Chisq(3), maxZgu-maxZunpaired)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    ret = string(ret, filename, " & Muse \$\\lambda_{\\text{GU}}=1\$ & ", @sprintf("%.2f", maxZgu), " & ",delta," & 13 & ", pvalstr, " & ", @sprintf("%.2f", maxZguparams[1]), " & ", @sprintf("%.2f", maxZguparams[2]), " & ", @sprintf("%.2f", maxZguparams[3]), "\\\\*\n")
    delta = @sprintf("%.2f", maxZ-maxZgu)
    pval = 1.0 - cdf(Chisq(1), maxZ-maxZgu)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    #H = @sprintf("%.1f", maxZentropy["H"])
    #H = @sprintf("%.1f", maxZentropy["H"])

    ret = string(ret, " & Muse unconstrained & ", @sprintf("%.2f", maxZ), " & ",delta," & 14 & ", pvalstr, " & ", @sprintf("%.2f", maxparams[1]), " & ", @sprintf("%.2f", maxparams[2]), " & ", @sprintf("%.2f", maxparams[3]), "\\tabularnewline\n")


    #\rowcolor{black!20} Unpaired & -15618.34 & - & 10 & - & 1.00 & 1.00 & 1.00 & 1.53 & 2.46 & 1.87 &  1.28 & 3.50 &  1.74\tabularnewline
    #Muse $\lambda_{\text{GU}}=1$ & -15423.56 & 194.78 & 13 & *** & 7.81 & 6.36 & 1.00 & 1.92 & 2.92 & 2.05 & 1.70 & 4.43 & 2.01\tabularnewline
    #\rowcolor{black!20} Muse unconstrained & -15419.26 & 4.30 & 14 & * & 10.63 & 8.88 & 2.16 & 2.01 & 2.85 & 2.13 & 1.73 & 4.28 & 2.09\tabularnewline
    #println(maxZ,"\t",maxparams)
    ret = string(ret, "\\midrule\n")
  end
  println(ret)
end

function pvaluehelper(alternativell, nullll, n1, n2, csvformat=false)
  delta = @sprintf("%.2f", alternativell-nullll)
  pval = 1.0 - cdf(Chisq(n1-n2), 2.0*(alternativell-nullll))
  pvalstr = "n.s."
  if pval < 0.05
    pvalstr = "*"
    #pvalstr = "p < 0.05"
  end
  if  pval < 0.005
    pvalstr = "**"
    #pvalstr = "p < 0.005"
  end
  if pval < 0.0005
    pvalstr = "***"
    #pvalstr = "p < 0.0005"
  end

  #pvalstr = string(@sprintf("%0.1e", pval), " (",pvalstr,")")
  if csvformat
      return string(delta,",", pvalstr)
  else
      return string(delta," & ", pvalstr)
  end
end

function printaictables2(files)
  ret = ""
  for maxfile in files
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]


    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    #maxZentropy = readentropy(string(maxfile, ".entropy"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))

    """
    delta = @sprintf("%.2f", maxZgu-maxZunpaired)
    pval = 1.0 - cdf(Chisq(3), maxZgu-maxZunpaired)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    """

    #=
    ret = string(ret, " & 1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    ret = string(ret, " & 2. Muse \$\\lambda_{\\text{GC}}=\\lambda_{\\text{AU}}, \\lambda_{\\text{GU}}=1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[1]), " & ", @sprintf("%.2f", maxZgcauguparams[2]), " & ", @sprintf("%.2f", maxZgcauguparams[3]), "\\\\*\n")
    ret = string(ret, filename, " & 3. Muse \$\\lambda_{\\text{GC}}=\\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,1}\$ = ",pvaluehelper(maxZgcau, maxZunpaired, 13, 10), " & ", @sprintf("%.2f", maxZgcauparams[1]), " & ", @sprintf("%.2f", maxZgcauparams[2]), " & ", @sprintf("%.2f", maxZgcauparams[3]), "\\\\*\n")
    ret = string(ret,  " & & & & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", " & ", " & ","\\\\*\n")
    ret = string(ret,  " & 4. Muse \$\\lambda_{\\text{GU}}=1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,1}\$ = ",pvaluehelper(maxZgu, maxZunpaired, 13, 10), " & ", @sprintf("%.2f", maxZguparams[1]), " & ", @sprintf("%.2f", maxZguparams[2]), " & ", @sprintf("%.2f", maxZguparams[3]), "\\\\*\n")
    ret = string(ret,  " & & & & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", " & ", " & ","\\\\*\n")
    =#


    ret = string(ret, " & 1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    ret = string(ret, " & 2. Muse \$\\lambda_{\\text{GC}}=\\lambda_{\\text{AU}}, \\lambda_{\\text{GU}}=1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[10]), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[19]), "\\\\*\n")
    ret = string(ret, filename, " & 3. Muse \$\\lambda_{\\text{GC}}=\\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,1}\$ = ",pvaluehelper(maxZgcau, maxZunpaired, 13, 10), " & ", @sprintf("%.2f", maxZgcauparams[10]), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[19]), "\\\\*\n")
    ret = string(ret,  " & & & & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", " & ", " & ","\\\\*\n")
    ret = string(ret,  " & 4. Muse \$\\lambda_{\\text{GU}}=1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,1}\$ = ",pvaluehelper(maxZgu, maxZunpaired, 13, 10), " & ", @sprintf("%.2f", maxZguparams[10]), " & ", @sprintf("%.2f", maxZguparams[17]), " & ", @sprintf("%.2f", maxZguparams[19]), "\\\\*\n")
    ret = string(ret,  " & & & & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", " & ", " & ","\\\\*\n")

    delta = @sprintf("%.2f", maxZ-maxZgu)
    pval = 1.0 - cdf(Chisq(1), maxZ-maxZgu)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    #H = @sprintf("%.1f", maxZentropy["H"])
    #H = @sprintf("%.1f", maxZentropy["H"])

    #ret = string(ret, " & 5. Muse unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,2}\$ = ",pvaluehelper(maxZ, maxZgcaugu, 14, 12), " & ", @sprintf("%.2f", maxparams[1]), " & ", @sprintf("%.2f", maxparams[2]), " & ", @sprintf("%.2f", maxparams[3]), "\\tabularnewline\n")
    ret = string(ret, " & 5. Muse unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,2}\$ = ",pvaluehelper(maxZ, maxZgcaugu, 14, 12), " & ", @sprintf("%.2f", maxparams[10]), " & ", @sprintf("%.2f", maxparams[17]), " & ", @sprintf("%.2f", maxparams[19]), "\\tabularnewline\n")
    ret = string(ret,  " & & & & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", " & ", " & ","\\\\*\n")
    ret = string(ret,  " & & & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\\\*\n")

    #\rowcolor{black!20} Unpaired & -15618.34 & - & 10 & - & 1.00 & 1.00 & 1.00 & 1.53 & 2.46 & 1.87 &  1.28 & 3.50 &  1.74\tabularnewline
    #Muse $\lambda_{\text{GU}}=1$ & -15423.56 & 194.78 & 13 & *** & 7.81 & 6.36 & 1.00 & 1.92 & 2.92 & 2.05 & 1.70 & 4.43 & 2.01\tabularnewline
    #\rowcolor{black!20} Muse unconstrained & -15419.26 & 4.30 & 14 & * & 10.63 & 8.88 & 2.16 & 2.01 & 2.85 & 2.13 & 1.73 & 4.28 & 2.09\tabularnewline
    #println(maxZ,"\t",maxparams)
    ret = string(ret, "\\midrule\n")
  end
  println(ret)
end

function printaictables(files)
  ret = ""
  for maxfile in files
    println(maxfile)
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]


    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    #maxZentropy = readentropy(string(maxfile, ".entropy"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))
    #ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", 1.00), " & ", @sprintf("%.2f", 1.00), " & ", @sprintf("%.2f", 1.00), "\\\\*\n")
    """
    delta = @sprintf("%.2f", maxZgu-maxZunpaired)
    pval = 1.0 - cdf(Chisq(3), maxZgu-maxZunpaired)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    """
    #ret = string(ret, "2. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,} \\lambda_{\\text{AU}}, \\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[1]), " & ", @sprintf("%.2f", maxZgcauguparams[2]), " & ", @sprintf("%.2f", maxZgcauguparams[3]), "\\\\*\n")
    #ret = string(ret, "3. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,}  \\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZgcauparams[1]), " & ", @sprintf("%.2f", maxZgcauparams[2]), " & ", @sprintf("%.2f", maxZgcauparams[3]), "\\\\*\n")
    #ret = string(ret, "4. \$\\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZguparams[1]), " & ", @sprintf("%.2f", maxZguparams[2]), " & ", @sprintf("%.2f", maxZguparams[3]), "\\\\*\n")
    ret = string(ret, "2. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,} \\lambda_{\\text{AU}}, \\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", 1.00), "\\\\*\n")
    ret = string(ret, "3. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,}  \\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[19]), "\\\\*\n")
    ret = string(ret, "4. \$\\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZguparams[10]), " & ", @sprintf("%.2f", maxZguparams[17]), " & ", @sprintf("%.2f", 1.00), "\\\\*\n")
    delta = @sprintf("%.2f", maxZ-maxZgu)
    pval = 1.0 - cdf(Chisq(1), maxZ-maxZgu)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    #H = @sprintf("%.1f", maxZentropy["H"])
    #H = @sprintf("%.1f", maxZentropy["H"])


    #ret = string(ret, "5. Unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", @sprintf("%.2f", maxparams[1]), " & ", @sprintf("%.2f", maxparams[2]), " & ", @sprintf("%.2f", maxparams[3]), "\\\\*\n")
    #ret = string(ret,  "& & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\tabularnewline\n")
    ret = string(ret, "5. Unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", @sprintf("%.2f", maxparams[10]), " & ", @sprintf("%.2f", maxparams[17]), " & ", @sprintf("%.2f", maxparams[19]), "\\\\*\n")
    ret = string(ret,  "& & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\tabularnewline\n")

    #\rowcolor{black!20} Unpaired & -15618.34 & - & 10 & - & 1.00 & 1.00 & 1.00 & 1.53 & 2.46 & 1.87 &  1.28 & 3.50 &  1.74\tabularnewline
    #Muse $\lambda_{\text{GU}}=1$ & -15423.56 & 194.78 & 13 & *** & 7.81 & 6.36 & 1.00 & 1.92 & 2.92 & 2.05 & 1.70 & 4.43 & 2.01\tabularnewline
    #\rowcolor{black!20} Muse unconstrained & -15419.26 & 4.30 & 14 & * & 10.63 & 8.88 & 2.16 & 2.01 & 2.85 & 2.13 & 1.73 & 4.28 & 2.09\tabularnewline
    #println(maxZ,"\t",maxparams)
    #ret = string(ret, "\\midrule\n")
  end
  println(ret)
  return ret
end

function printaictablescsv(files)
  ret = ""
  for maxfile in files
    println(maxfile)
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]

    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))
    ret = "Model,Params,Maximum Log Likelihood,Delta,p-value,lambdaGC,lambdaAU,lambdaGU\n"
    ret = string(ret, "\"1. Unpaired\", 10,", @sprintf("%.2f", maxZunpaired), ",","",",", "-,", @sprintf("%.2f", 1.00), ",", @sprintf("%.2f", 1.00), ",", @sprintf("%.2f", 1.00), "\n")
    ret = string(ret, "\"2. lambdaGC := lambdaAU, lambdaGU := 1\",12,", @sprintf("%.2f", maxZgcaugu), ",M2-M1 = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10, true), ",", @sprintf("%.2f", maxZgcauguparams[17]), ",", @sprintf("%.2f", maxZgcauguparams[17]), ",", @sprintf("%.2f", 1.00), "\n")
    ret = string(ret, "\"3. lambdaGC := lambdaAU\",13,", @sprintf("%.2f", maxZgcau), ",M3-M2 = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12, true), ",", @sprintf("%.2f", maxZgcauparams[17]), ",", @sprintf("%.2f", maxZgcauparams[17]), ",", @sprintf("%.2f", maxZgcauparams[19]), "\n")
    ret = string(ret, "\"4. lambdaGU := 1\",13,", @sprintf("%.2f", maxZgu), ",M4-M2 = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12, true), ",", @sprintf("%.2f", maxZguparams[10]), ",", @sprintf("%.2f", maxZguparams[17]), ",", @sprintf("%.2f", 1.00), "\n")
    ret = string(ret, "\"5. Unconstrained\",14,", @sprintf("%.2f", maxZ), ",M5-M3 = ",pvaluehelper(maxZ, maxZgcau, 14, 13, true), ",", @sprintf("%.2f", maxparams[10]), ",", @sprintf("%.2f", maxparams[17]), ",", @sprintf("%.2f", maxparams[19]), "\n")
    ret = string(ret,  ",,,M5-M4 = ",pvaluehelper(maxZ, maxZgu, 14, 13, true), ",", ",", ",","\n")
    ret = string(ret, "\"*p < 0.05; **p < 0.005; ***p < 0.0005; n.s. = not significant\"\n")
  end
  return ret
end

function printlargeaictable(files)
  ret = ""
  for maxfile in files
    println(maxfile)
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]


    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    #maxZentropy = readentropy(string(maxfile, ".entropy"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))
    #ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[10]), " & ", @sprintf("%.2f", maxunpairedparams[17]), " & ", @sprintf("%.2f", maxunpairedparams[19]), "\\\\*\n")

    ret = string(ret, "2. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,} \\lambda_{\\text{AU}}, \\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[19]), "\\\\*\n")
    ret = string(ret, "3. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,}  \\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[19]), "\\\\*\n")
    ret = string(ret, "4. \$\\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZguparams[10]), " & ", @sprintf("%.2f", maxZguparams[17]), " & ", @sprintf("%.2f", maxZguparams[19]), "\\\\*\n")
    delta = @sprintf("%.2f", maxZ-maxZgu)
    pval = 1.0 - cdf(Chisq(1), maxZ-maxZgu)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end

    ret = string(ret, "5. Unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", @sprintf("%.2f", maxparams[10]), " & ", @sprintf("%.2f", maxparams[17]), " & ", @sprintf("%.2f", maxparams[19]), "\\\\*\n")
    ret = string(ret,  "& & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\tabularnewline\n")

    ret = string(ret, "\\midrule\n")
  end
  println(ret)
end

function printaictables3(files)
  ret = ""
  for maxfile in files
    println(maxfile)
    m = match(r"^(.*/)([^/]*)$", maxfile)
    filename = m[2]


    maxZ, maxparams = readmaxparams(string(maxfile, ".max"))
    #maxZentropy = readentropy(string(maxfile, ".entropy"))
    maxZgu, maxZguparams = readmaxparams(string(maxfile, ".fixgu.max"))
    maxZgcau, maxZgcauparams = readmaxparams(string(maxfile, ".fixgcau.max"))
    maxZgcaugu, maxZgcauguparams = readmaxparams(string(maxfile, ".fixgcaugu.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(maxfile, ".max.unpaired"))
    #ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[1]), " & ", @sprintf("%.2f", maxunpairedparams[2]), " & ", @sprintf("%.2f", maxunpairedparams[3]), "\\\\*\n")
    ret = string(ret, "1. Unpaired &  10 & ", @sprintf("%.2f", maxZunpaired), " & ",""," &", " - & ", @sprintf("%.2f", maxunpairedparams[10]), " & ", @sprintf("%.2f", maxunpairedparams[17]), " & ", @sprintf("%.2f", maxunpairedparams[19]), "\\\\*\n")
    """
    delta = @sprintf("%.2f", maxZgu-maxZunpaired)
    pval = 1.0 - cdf(Chisq(3), maxZgu-maxZunpaired)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    """
    #ret = string(ret, "2. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,} \\lambda_{\\text{AU}}, \\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[1]), " & ", @sprintf("%.2f", maxZgcauguparams[2]), " & ", @sprintf("%.2f", maxZgcauguparams[3]), "\\\\*\n")
    #ret = string(ret, "3. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,}  \\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZgcauparams[1]), " & ", @sprintf("%.2f", maxZgcauparams[2]), " & ", @sprintf("%.2f", maxZgcauparams[3]), "\\\\*\n")
    #ret = string(ret, "4. \$\\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZguparams[1]), " & ", @sprintf("%.2f", maxZguparams[2]), " & ", @sprintf("%.2f", maxZguparams[3]), "\\\\*\n")
    ret = string(ret, "2. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,} \\lambda_{\\text{AU}}, \\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 12 & ", @sprintf("%.2f", maxZgcaugu), " & \$\\Delta_{2,1}\$ = ", pvaluehelper(maxZgcaugu, maxZunpaired, 12, 10), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[17]), " & ", @sprintf("%.2f", maxZgcauguparams[19]), "\\\\*\n")
    ret = string(ret, "3. \$\\lambda_{\\text{GC}} {\\,\\coloneqq\\,}  \\lambda_{\\text{AU}}\$ &  13 & ", @sprintf("%.2f", maxZgcau), " & \$\\Delta_{3,2}\$ = ",pvaluehelper(maxZgcau, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[17]), " & ", @sprintf("%.2f", maxZgcauparams[19]), "\\\\*\n")
    ret = string(ret, "4. \$\\lambda_{\\text{GU}} {\\,\\coloneqq\\,}  1\$ & 13 &", @sprintf("%.2f", maxZgu), " & \$\\Delta_{4,2}\$ = ",pvaluehelper(maxZgu, maxZgcaugu, 13, 12), " & ", @sprintf("%.2f", maxZguparams[10]), " & ", @sprintf("%.2f", maxZguparams[17]), " & ", @sprintf("%.2f", maxZguparams[19]), "\\\\*\n")
    delta = @sprintf("%.2f", maxZ-maxZgu)
    pval = 1.0 - cdf(Chisq(1), maxZ-maxZgu)
    pvalstr = "n.s."
    if pval < 0.01
      pvalstr = "*"
    end
    if  pval < 0.001
      pvalstr = "**"
    end
    if pval < 0.0001
      pvalstr = "***"
    end
    #H = @sprintf("%.1f", maxZentropy["H"])
    #H = @sprintf("%.1f", maxZentropy["H"])


    #ret = string(ret, "5. Unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", @sprintf("%.2f", maxparams[1]), " & ", @sprintf("%.2f", maxparams[2]), " & ", @sprintf("%.2f", maxparams[3]), "\\\\*\n")
    #ret = string(ret,  "& & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\tabularnewline\n")
    ret = string(ret, "5. Unconstrained & 14 & ", @sprintf("%.2f", maxZ), " & \$\\Delta_{5,3}\$ = ",pvaluehelper(maxZ, maxZgcau, 14, 13), " & ", @sprintf("%.2f", maxparams[10]), " & ", @sprintf("%.2f", maxparams[17]), " & ", @sprintf("%.2f", maxparams[19]), "\\\\*\n")
    ret = string(ret,  "& & & \$\\Delta_{5,4}\$ = ",pvaluehelper(maxZ, maxZgu, 14, 13), " & ", " & ", " & ","\\tabularnewline\n")

    #\rowcolor{black!20} Unpaired & -15618.34 & - & 10 & - & 1.00 & 1.00 & 1.00 & 1.53 & 2.46 & 1.87 &  1.28 & 3.50 &  1.74\tabularnewline
    #Muse $\lambda_{\text{GU}}=1$ & -15423.56 & 194.78 & 13 & *** & 7.81 & 6.36 & 1.00 & 1.92 & 2.92 & 2.05 & 1.70 & 4.43 & 2.01\tabularnewline
    #\rowcolor{black!20} Muse unconstrained & -15419.26 & 4.30 & 14 & * & 10.63 & 8.88 & 2.16 & 2.01 & 2.85 & 2.13 & 1.73 & 4.28 & 2.09\tabularnewline
    #println(maxZ,"\t",maxparams)
    ret = string(ret, "\\midrule\n")
  end
  println(ret)
end

function printcomputationtables(files)
  iter = 1
  ret = ""
  for countfile in files
    fin = open(string(countfile,".calculations"), "r")
    line1 = ""
    line2 = ""
    i = 1
    for line in readlines(fin)
      if i == 1
        line1 = strip(line)
      elseif i == 2
        line2 = strip(line)
      else
        break
      end
      i += 1
    end
    close(fin)
    spl1 = split(strip(line1,['(',')']),",")
    spl2 = split(line2)
    without = parse(Int, string(spl1[4]))
    with = parse(Int, string(spl1[3]))
    m = match(r"^(.*/)([^/]*)$", countfile)
    if iter % 2 == 1
      ret = string(ret, "\\rowcolor{black!20} ")
    end
    ret = string(ret,m[2], " & ", spl2[1], " & ", spl2[2], " & ", format(without, commas=true), " & ", format(with, commas=true), " & ", @sprintf("%.1f", without/with), "\\tabularnewline\n")
    iter += 1
  end
  println(ret)
end

function printcomputationtablessorted(files)
  iter = 1
  ret = ""
  lines = []
  for countfile in files
    fin = open(string(countfile,".calculations"), "r")
    line1 = ""
    line2 = ""
    i = 1
    for line in readlines(fin)
      if i == 1
        line1 = strip(line)
      elseif i == 2
        line2 = strip(line)
      else
        break
      end
      i += 1
    end
    close(fin)
    spl1 = split(strip(line1,['(',')']),",")
    spl2 = split(line2)
    without = parse(Int, string(spl1[4]))
    with = parse(Int, string(spl1[3]))
    m = match(r"^(.*/)([^/]*)$", countfile)
    ret = ""
    #=
    if iter % 2 == 1
      ret = string(ret, "\\rowcolor{black!20} ")
    end=#
    ret = string(ret,m[2], " & ", spl2[1], " & ", spl2[2], " & ", format(without, commas=true), " & ", format(with, commas=true), " & ", @sprintf("%.1f", without/with), "\\tabularnewline\n")
    push!(lines, (-without, ret))
    iter += 1
  end
  sort!(lines)
  i = 1
  for line in lines
    if i % 2 == 1
      print("\\rowcolor{black!20} ")
    end
    print(line[2])
    i += 1
  end

  #println(ret)
end

function readmaxparams(maxfile)
  jsondict = JSON.parse(open(maxfile))
  return (jsondict["Z"], convert(Array{Float64,1}, jsondict["maxparams"]))
end

function readentropy(entropyfile)
  return JSON.parse(open(entropyfile))
end

function printentropytable(files)
  ret = ""
  iter = 1
  for entropyfile in files
    m = match(r"^(.*/)([^/]*)$", entropyfile)
    filename = m[2]
    maxZ, maxparams = readmaxparams(string(entropyfile, ".max"))
    entropydict = readentropy(string(entropyfile, ".entropy"))
    maxZ_s1, maxparams_s1 = readmaxparams(string(entropyfile, ".shuffle1.max"))
    #entropypriordict = readentropy(string(entropyfile, ".entropyprior"))
    entropydict_s1 = readentropy(string(entropyfile, ".shuffle1.entropy"))
    entropydict_s2 = readentropy(string(entropyfile, ".shuffle2.entropy"))
    maxZ_s2, maxparams_s2 = readmaxparams(string(entropyfile, ".shuffle2.max"))
    maxZunpaired, maxunpairedparams = readmaxparams(string(entropyfile, ".max.unpaired"))
    #println(entropydict)
    #println(entropydict_s1)
    if iter % 2 == 1
      ret = string(ret, "\\rowcolor{black!20} ")
    end
    #ret = string(ret, filename, " & ",  @sprintf("%.1f", maxZ), " & ", entropydict["length"], " & ", @sprintf("%.1f", entropydict["H"]), " & ", @sprintf("%.1f", entropydict["Hmax"]), " & ", @sprintf("%.1f", entropypriordict["H"]), " & ", @sprintf("%.3f", entropydict["percentage"]), " & ",  @sprintf("%.1f", maxZ_s1), " & ", @sprintf("%.1f", entropydict_s1["H"]), " & ", @sprintf("%.3f", entropydict_s1["percentage"]), " & ", @sprintf("%.1f", maxZ_s2), " & ", @sprintf("%.1f", entropydict_s2["H"]), " & ", @sprintf("%.3f", entropydict_s2["percentage"]), "\\tabularnewline\n")
    #ret = string(ret, filename, " & ",  @sprintf("%.1f", maxZ), " & ", entropydict["length"], " & ", @sprintf("%.1f", entropydict["H"]), " & ", @sprintf("%.1f", entropydict["Hmax"]), " & ", @sprintf("%.1f", entropypriordict["H"]), " & ", @sprintf("%.3f", entropydict["percentage"]), " & ",  @sprintf("%.1f", maxZ_s1-maxZ), " & ", @sprintf("%.1f", entropydict_s1["H"]-entropydict["H"]), " & ", @sprintf("%.3f", entropydict_s1["percentage"]), " & ", @sprintf("%.1f", maxZ_s2-maxZ), " & ", @sprintf("%.1f", entropydict_s2["H"]-entropydict["H"]), " & ", @sprintf("%.3f", entropydict_s2["percentage"]), "\\tabularnewline\n")
    #ret = string(ret, filename, " & ",  @sprintf("%.1f", maxZ), " & ", entropydict["length"], " & ", @sprintf("%.1f", entropydict["H"]), " & ", @sprintf("%.1f", entropydict["Hmax"]), " & ", @sprintf("%.3f", entropydict["percentage"]), " & ",  @sprintf("%.1f", maxZ_s1-maxZ), " & ", @sprintf("%.1f", entropydict_s1["H"]-entropydict["H"]), " & ", @sprintf("%.1f", (entropydict_s1["H"]-entropydict["H"])*100.0/entropydict["Hmax"]), "\\%", " & ", @sprintf("%.1f", maxZ_s2-maxZ), " & ", @sprintf("%.1f", entropydict_s2["H"]-entropydict["H"]), " & ", @sprintf("%.1f", (entropydict_s2["H"]-entropydict["H"])*100.0/entropydict["Hmax"]),  "\\%", "\\tabularnewline\n")

    deltaH1 = @sprintf("%.1f", entropydict_s1["H"]-entropydict["H"])
    deltaH1perc = @sprintf("%.3f", (entropydict_s1["H"]-entropydict["H"])/entropydict["Hmax"])
    if entropydict_s1["H"]-entropydict["H"] >= 0.0
      deltaH1 = string("+", @sprintf("%.1f", entropydict_s1["H"]-entropydict["H"]))
      deltaH1perc = string("+", @sprintf("%.3f", (entropydict_s1["H"]-entropydict["H"])/entropydict["Hmax"]))
    end

    deltaH2 = @sprintf("%.1f", entropydict_s2["H"]-entropydict["H"])
    deltaH2perc = @sprintf("%.3f", (entropydict_s2["H"]-entropydict["H"])/entropydict["Hmax"])
    if entropydict_s2["H"]-entropydict["H"] >= 0.0
      deltaH2 = string("+", @sprintf("%.1f", entropydict_s2["H"]-entropydict["H"]))
      deltaH2perc = string("+", @sprintf("%.3f", (entropydict_s2["H"]-entropydict["H"])/entropydict["Hmax"]))
    end


    ret = string(ret, filename, " & ",  @sprintf("%.1f", maxZ), " & ", @sprintf("%.1f", maxZunpaired-maxZ), " & ", entropydict["length"], " & ", @sprintf("%.1f", entropydict["H"]), " & ", @sprintf("%.1f", entropydict["Hmax"]), " & ", @sprintf("%.3f", entropydict["percentage"]), " & ",  @sprintf("%.1f", maxZ_s1-maxZ), " & ", deltaH1, " & ", deltaH1perc, "", " & ", @sprintf("%.1f", maxZ_s2-maxZ), " & ", deltaH2, " & ", deltaH2perc,  "", "\\tabularnewline\n")
    iter += 1
  end
  println(ret)
end

function getcovmatfromfile(mcmclogfile)
  fin = open(mcmclogfile,"r")
  params = Dict{AbstractString,Array{Float64,1}}()
  header = true
  keys = AbstractString[]
  for line in readlines(fin)
    spl = split(strip(line))
    if header
      keys = spl
      for key in keys
        params[key] = Float64[]
      end
      header = false
    else
      for i=1:length(spl)
        p = get(params, keys[i], Float64[])
        push!(p, parse(Float64,spl[i]))
        params[keys[i]] = p
      end
    end
  end
  close(fin)
  len = -1
  for key in keys
    len = length(params[key])
    if len > 1
      data = params[key][div(len,2):end]
      println(key,"\t",length(data),"\t", mean(data), "\t", std(data))
    end
  end

  tuningvectors = Array{Float64,1}[]
  startindex = div(len,2)
  len2 = length(params["lambdaGC"])
  println("LEN",len2)
  for i=startindex:len2
    tvec = Float64[params["lambdaGC"][i], params["lambdaAT"][i], params["lambdaGT"][i], params["qAC"][i], params["qAG"][i], params["qAT"][i], params["qCG"][i], params["qCT"][i], params["qGT"][i], params["lambdaGammaShape"][i], params["siteGammaShape"][i], params["lambdazeroweight"][i]]
    push!(tuningvectors,tvec)
  end


  mat = hcat(tuningvectors[1],tuningvectors[2])
  for z=3:length(tuningvectors)
    mat = hcat(mat,tuningvectors[z])
  end
  mat = mat'

  #println(mat)
  #covmat = (2.38*2.38*cov(mat))/d2

  #=
  len = length(params[keys[1]])
  println("lambdaGC > lambdaAT: ", mean([lambdaGC > lambdaAT ? 1.0 : 0.0 for (lambdaAT,lambdaGC) in zip(params["lambdaAT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))
  println("lambdaAT > lambdaGT: ", mean([lambdaAT > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaAT) in zip(params["lambdaGT"][div(len,2):end],params["lambdaAT"][div(len,2):end])]))
  println("lambdaGC > lambdaGT: ", mean([lambdaGC > lambdaGT ? 1.0 : 0.0 for (lambdaGT,lambdaGC) in zip(params["lambdaGT"][div(len,2):end],params["lambdaGC"][div(len,2):end])]))

  startindex = div(len,2)
  len2 = length(params["lambdaGC"][startindex:end])
  data = zeros(len2,3)
  for i=1:len2
    data[i,1] = params["lambdaGC"][startindex+i-1]
    data[i,2] = params["lambdaAT"][startindex+i-1]
    data[i,3] = params["lambdaGT"][startindex+i-1]
  end
  A = full(chol(cov(data)))
  return keys,params,A=#
  return cov(mat)
end

function saveungappedlength(outputprefix)
  fastafile = string(outputprefix,".fas.norm")
  nucmapping = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'U' => 4)
  sequences = AbstractString[]
  names = AbstractString[]
  seqnametoindex = Dict{AbstractString,Int}()
  len = 0
  seqindex = 1
  FastaIO.FastaReader(fastafile) do fr
   for (desc, seq) in fr
     len = length(seq)
     push!(names,desc)
     push!(sequences, seq)
     seqnametoindex[desc] = seqindex
     seqindex += 1
   end
  end
  numseqs = length(sequences)
  numcols = len
  data = zeros(Float64,numseqs,numcols,4)
  obsfreqs = zeros(Float64,4)
  gapfrequency = zeros(Float64,numcols)
  for s=1:numseqs
    seq = sequences[s]
    for j=1:numcols
      nuc = get(nucmapping,seq[j],0)
      if nuc > 0
        data[s,j,nuc] = 1.0
        obsfreqs[nuc] += 1.0
      else
        data[s,j,:] = 1.0
        gapfrequency[j] += 1.0
      end
    end
  end
  gapfrequency /= numseqs

  cutoff = 0.5
  index = 1
  for i=1:numcols
    if gapfrequency[i] < cutoff
      index += 1
    end
  end
  newlen = index - 1

  println("NEWLENGTH ", newlen,"\t", numcols)
  fout = open(string(outputprefix,".length"),"w")
  write(fout, string(newlen,"\n"))
  write(fout, string(numseqs,"\n"))
  close(fout)
end

function lambdaordering(outputprefix)
    params,len = summarizemcmc(string(outputprefix,".B1.0.M1.0.mcmc.log"))
    GCvec = params["lambdaGC"][div(len,2):end]
    AUvec = params["lambdaAT"][div(len,2):end]
    GUvec = params["lambdaGT"][div(len,2):end]

    counts = zeros(Float64, 6)
    for (GC,AU,GU) in zip(GCvec,AUvec,GUvec)
      if GC >= AU >= GU
        counts[1] += 1.0
      elseif AU >= GC >= GU
        counts[2] += 1.0
      elseif GC >= GU >= AU
        counts[3] += 1.0
      elseif AU >= GU >= GC
        counts[4] += 1.0
      elseif GU >= GC >= AU
        counts[5] += 1.0
      else # GU >= AU >= GC
        counts[6] += 1.0
      end
    end

    return counts/sum(counts), median(GCvec), median(AUvec), median(GUvec)
end

function main()
    plotstructurebenchmarks()
    exit()

    #files = ["/media/michael/Sandisk500GB/data/msv/msv","/media/michael/Sandisk500GB/data/tylcv/tylcv","/media/michael/Sandisk500GB/data4/fmdv2/fmdv","/media/michael/Sandisk500GB/data4/hpv1/hpv1","/media/michael/Sandisk500GB/data4/tobamovirus/tobamovirus","/media/michael/Sandisk500GB/data4/rhinovirus_a/rhinovirus_a","/media/michael/Sandisk500GB/data4/hepa/hepa","/media/michael/Sandisk500GB/data/RNaseMRP_mcmc/RNaseMRP", "/media/michael/Sandisk500GB/data/ires_mcmc6/ires", "/media/michael/Sandisk500GB/data/RF00001/RF00001","/media/michael/Sandisk500GB/data/RF00002/RF00002","/media/michael/Sandisk500GB/data/RF00003/RF00003", "/media/michael/Sandisk500GB/data/RF00004/RF00004", "/media/michael/Sandisk500GB/data/RF00010/RF00010", "/media/michael/Sandisk500GB/data/RF00011/RF00011", "/media/michael/Sandisk500GB/data/RF00012/RF00012", "/media/michael/Sandisk500GB/data/RF00020/RF00020", "/media/michael/Sandisk500GB/data/RF00026/RF00026", "/media/michael/Sandisk500GB/data/RF00100/RF00100", "/media/michael/Sandisk500GB/data/RF00174/RF00174", "/media/michael/Sandisk500GB/data/RF00379/RF00379", "/media/michael/Sandisk500GB/data/RF00380/RF00380", "/media/michael/Sandisk500GB/data/RF01846/RF01846", "/media/michael/Sandisk500GB/data/RF01854/RF01854", "/media/michael/Sandisk500GB/data/RF02001/RF02001", "/media/michael/Sandisk500GB/data/RF02540/RF02540", "/media/michael/Sandisk500GB/data/RF02541/RF02541", "/media/michael/Sandisk500GB/data/RF02542/RF02542", "/media/michael/Sandisk500GB/data/RF02543/RF02543"  , "/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv", "/media/michael/Sandisk500GB/data/wdf_mcmc/wdf","/media/michael/Sandisk500GB/data/bocavirus_mcmc/bocavirus"]
    #printcomputationtables(files)
    #printcomputationtablessorted(files)
    #exit()

    #println(getcovmatfromfile("/media/michael/Sandisk500GB/data/wdf_mcmc/wdf.B1.0.M1.0.mcmc.log"))
    #exit()

    #files = ["/media/michael/Sandisk500GB/data/RF00010_mcmc/RF00010", "/media/michael/Sandisk500GB/data/wdf_mcmc/wdf", "/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv","/media/michael/Sandisk500GB/data/RF00001/RF00001","/media/michael/Sandisk500GB/data/RF00002/RF00002","/media/michael/Sandisk500GB/data/RF00003/RF00003"]
    #files = ["/media/michael/Sandisk500GB/data/RF01846/RF01846"]

    files = ["/media/michael/Sandisk500GB/data4/beakandfeatherdisease/beakandfeatherdisease","/media/michael/Sandisk500GB/data4/fmdv2/fmdv","/media/michael/Sandisk500GB/data4/hpv1/hpv1","/media/michael/Sandisk500GB/data4/tobamovirus/tobamovirus","/media/michael/Sandisk500GB/data4/rhinovirus_a/rhinovirus_a","/media/michael/Sandisk500GB/data4/hepa/hepa","/media/michael/Sandisk500GB/data/wdf_mcmc/wdf","/media/michael/Sandisk500GB/data/bocavirus_mcmc/bocavirus","/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv","/media/michael/Sandisk500GB/data/tylcv/tylcv","/media/michael/Sandisk500GB/data/msv/msv"]
    files = files[1:1]
    #printaictables(files)
    printlargeaictable(files)
    exit()
    #=
    files = ["/media/michael/Sandisk500GB/data/RF00002_mcmc/RF00002", "/media/michael/Sandisk500GB/data/RF00002_cpu/RF00002","/media/michael/Sandisk500GB/data/RF00003_mcmc/RF00003", "/media/michael/Sandisk500GB/data/RF00003_cpu/RF00003", "/media/michael/Sandisk500GB/data/RF00010_mcmc/RF00010_c2", "/media/michael/Sandisk500GB/data/RF00010_cpu/RF00010","/media/michael/Sandisk500GB/data/RF00011_mcmc/RF00011", "/media/michael/Sandisk500GB/data/RF00011_cpu/RF00011", "/media/michael/Sandisk500GB/data/RF00020_mcmc/RF00020", "/media/michael/Sandisk500GB/data/RF00020_cpu/RF00020", "/media/michael/Sandisk500GB/data/RF01846_mcmc/RF01846", "/media/michael/Sandisk500GB/data/RF01846_cpu/RF01846", "/media/michael/Sandisk500GB/data/RF00379_mcmc/RF00379", "/media/michael/Sandisk500GB/data/RF00379_cpu/RF00379", "/media/michael/Sandisk500GB/data/RF02542_mcmc/RF02542", "/media/michael/Sandisk500GB/data/RF02542_cpu/RF02542"]

    for file in files
      saveungappedlength(file)
    end

    exit()=#

    #files = ["/media/michael/Sandisk500GB/data/RF00001/RF00001","/media/michael/Sandisk500GB/data/RF00002/RF00002","/media/michael/Sandisk500GB/data/RF00003/RF00003", "/media/michael/Sandisk500GB/data/RF00004/RF00004", "/media/michael/Sandisk500GB/data/RF00026/RF00026", "/media/michael/Sandisk500GB/data/RF00100/RF00100", "/media/michael/Sandisk500GB/data/RF00174/RF00174", "/media/michael/Sandisk500GB/data/wdf/wdf"]
    #printtables(files)
    files = ["/media/michael/Sandisk500GB/data/wdf/wdf", "/media/michael/Sandisk500GB/data/msv/msv","/media/michael/Sandisk500GB/data/RF00001/RF00001", "/media/michael/Sandisk500GB/data/RF00002/RF00002", "/media/michael/Sandisk500GB/data/RF00003/RF00003", "/media/michael/Sandisk500GB/data/RF00004/RF00004","/media/michael/Sandisk500GB/data4/hepa/hepa", "/media/michael/Sandisk500GB/data4/rhinovirus_a/rhinovirus_a", "/media/michael/Sandisk500GB/data4/hpv1/hpv1"]
    #printentropytable(files)
    #exit()

    #outputprefix = "/media/michael/Sandisk500GB/data/RF00001_mcmc2/RF00001"
    #outputprefix = "/media/michael/Sandisk500GB/data/wdf_mcmc/wdf"
    #outputprefix = "/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv"

    #files = ["/media/michael/Sandisk500GB/data/RF00010_mcmc/RF00010","/media/michael/Sandisk500GB/data/RF00003_mcmc/RF00003", "/media/michael/Sandisk500GB/data/wdf_mcmc/wdf", "/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv","/media/michael/Sandisk500GB/data/RF00001_mcmc2/RF00001"]
    #files = ["/media/michael/Sandisk500GB/data/tylcv/tylcv"]
    #
    #"/media/michael/Sandisk500GB/data4/beakandfeatherdisease/beakandfeatherdisease"
    files = ["/media/michael/Sandisk500GB/data/wdf_mcmc/wdf","/media/michael/Sandisk500GB/data/bocavirus_mcmc/bocavirus","/media/michael/Sandisk500GB/data/beet_curly_mcmc/bctv","/media/michael/Sandisk500GB/data/tylcv/tylcv","/media/michael/Sandisk500GB/data/msv/msv", "/media/michael/Sandisk500GB/data/RF00001_mcmc2/RF00001","/media/michael/Sandisk500GB/data/RF00003_mcmc/RF00003", "/media/michael/Sandisk500GB/data/RF00010_mcmc/RF00010", "/media/michael/Sandisk500GB/data/RF00379_mcmc/RF00379", "/media/michael/Sandisk500GB/data/RF01846_mcmc/RF01846","/media/michael/Sandisk500GB/data4/hepa/hepa","/media/michael/Sandisk500GB/data4/rhinovirus_a/rhinovirus_a","/media/michael/Sandisk500GB/data4/tobamovirus/tobamovirus","/media/michael/Sandisk500GB/data4/hpv1/hpv1","/media/michael/Sandisk500GB/data4/fmdv2/fmdv"]

    data = []
    names = []
    medians = []
    for outputprefix in files
      println(outputprefix)
      m = match(r"^(.*/)([^/]*)$", outputprefix)
      push!(names,m[2])
      #outputdir = m[1]
      proportions, medianGC, medianAU, medianGU = lambdaordering(outputprefix)
      push!(data, proportions)
      println(proportions)
      #plotdistributions(outputprefix)
      #plotratios(outputprefix)
      push!(medians, [medianGC, medianAU, medianGU])
    end
    println(names)
    println(data)
    println(medians)
    #
end

function calculateGCcontent(sequences::Array{AbstractString,1})
  GC = 0.0
  AT = 0.0
  for seq in sequences
    for c in uppercase(seq)
      if c == 'G' || c == 'C'
        GC += 1.0
      elseif c == 'A' || c == 'T'  || c == 'U'
        AT += 1.0
      end
    end
  end
  return GC/(GC+AT)
end

function calculategapfrequency(sequences::Array{AbstractString,1})
  gaps = 0.0
  total = 0.0
  for seq in sequences
    for c in uppercase(seq)
      if c == '-'
        gaps += 1.0
      end
      total += 1.0
    end
  end
  return gaps/total
end

function calculatemediangapfrequency(sequences::Array{AbstractString,1})
  gapfreqs = Float64[]
  for col=1:length(sequences[1])
    gaps = 0.0
    total = 0.0
    for seq in sequences
      if seq[col] == '-'
        gaps += 1.0
      end
      total += 1.0
    end
    push!(gapfreqs, gaps/total)
  end
  return median(gapfreqs)
end

function calculateavgpairwiseidentity(sequences::Array{AbstractString,1})
  diversities = Float64[]
  for (i,s1) in enumerate(sequences)
    for (j,s2) in enumerate(sequences)
      if i > j
        matches = 0.0
        total = 0.0
        for (a,b) in zip(s1,s2)
          if a != '-' && b != '-'
            if a == b
              matches += 1.0
            end
            total += 1.0
          end
        end
        pairwisediversity = (matches/total)
        push!(diversities, pairwisediversity)
      end
    end
  end
  return mean(diversities)
end

function printalignmentsummarytable(alignmentfiles::Array{AbstractString,1})
  for alignment in alignmentfiles
    m = match(r"([^\.]+)\.?.*", basename(alignment))
    datasetname = m[1]
    
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(alignment) do fr
        for (desc, seq) in fr
           push!(names,desc)
           push!(sequences, seq)
        end
    end
    GCcontent = calculateGCcontent(sequences)
    gapfrequency = calculategapfrequency(sequences)
    mediancolumngapfrequency = calculatemediangapfrequency(sequences)
    avgpairwiseidentity = calculateavgpairwiseidentity(sequences)
    println(datasetname," & ", length(sequences)," & ", length(sequences[1]), " & ", @sprintf("%.1f", GCcontent*100.0),"%", " & ", @sprintf("%.1f", gapfrequency*100.0),"%", " & ", @sprintf("%.1f", mediancolumngapfrequency*100.0),"%", " & ", @sprintf("%.1f", avgpairwiseidentity*100.0), "%")
  end
end
#=
rfamdir = "../datasets/RFAM/3D/"
alignmentfiles = AbstractString[]
for f in filter(x -> match(r".+\.fas", x) != nothing, readdir(rfamdir))  
  push!(alignmentfiles, joinpath(rfamdir,f))
end
println("numfiles $(length(alignmentfiles))")

printalignmentsummarytable(alignmentfiles)
=#
#main()
