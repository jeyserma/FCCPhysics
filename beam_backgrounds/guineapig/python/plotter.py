import ROOT
import math

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def make_dynamic_legend(
    hists,
    labels,
    xmin,
    xmax,
    marker_size=0.0,
    text_size=0.035,
    entry_height=0.045,
    margin=0.02,
    top=0.90,
    min_bottom=0.15,
    legend_header=None,
    legend_ncols=1
):
    """
    Create a dynamically positioned ROOT legend.

    Parameters
    ----------
    hists : list[TH1]
        Histograms (already styled).
    labels : list[str]
        Legend labels.
    xmin, xmax : float
        X-axis range (used to decide left/right placement).
    marker_size : float
        If >0, use marker+line legend entries.
    text_size : float
        Legend text size.
    entry_height : float
        Vertical space per entry (NDC).
    margin : float
        Extra padding.
    top : float
        Top position (NDC).
    min_bottom : float
        Minimum bottom margin (avoid overlap with axis).

    Returns
    -------
    ROOT.TLegend
    """

    n_entries = len(hists)
    if legend_header:
        n_entries += 1
    if legend_ncols > 1:
        n_entries = math.ceil(1.0*n_entries/legend_ncols)+1

    # --- compute legend height ---
    leg_height = n_entries * entry_height + margin



    # --- determine where the "bulk" of the data is ---
    # Use first histogram as heuristic (fast & robust enough)
    h0 = hists[0]
    max_bin = h0.GetMaximumBin()
    x_peak = h0.GetXaxis().GetBinCenter(max_bin)

    x_mid = 0.5 * (xmin + xmax)

    # --- choose left or right ---
    if x_peak > x_mid:
        x1, x2 = 0.18, 0.45
    else:
        x1, x2 = 0.60, 0.88

    x1, x2 = 0.18, 0.85
    # --- vertical placement ---
    y2 = top
    y1 = y2 - leg_height

    # safety clamp
    if y1 < min_bottom:
        y1 = min_bottom
        y2 = y1 + leg_height

    # --- create legend ---
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(text_size)
    leg.SetNColumns(legend_ncols)
    if legend_header:
        leg.SetHeader(legend_header)

    # --- fill legend ---
    for h, lab in zip(hists, labels):
        opt = "l"
        if marker_size > 0:
            opt = "lp"
        leg.AddEntry(h, lab, opt)

    return leg

def plot_hists_1d(
    hists,
    labels,
    outname,
    x_title=None,
    y_title=None,
    x_range=None,          # tuple: (xmin, xmax)
    y_range=None,          # tuple: (ymin, ymax)
    logy=False,
    logx=False,
    colors=None,
    line_width=3,
    line_styles=None,
    marker_styles=None,
    marker_size=0.0,
    legend_pos=(0.62, 0.70, 0.88, 0.88),
    legend_header=None,
    legend_ncols=1,
    canvas_size=(800, 800),
    normalize=False,
    draw_option="hist",
    extra_text_left=None,
    extra_text_right=None,
    rebin=None,
    scale_factor=None,
    y_factor=1.35,
    print_bins=False,
):


    if len(hists) == 0:
        raise ValueError("No histograms provided.")
    if len(labels) > 0 and len(hists) != len(labels):
        raise ValueError("hists and labels must have the same length.")

    ROOT.gStyle.SetOptStat(0)

    # Default colors
    if colors is None:
        colors = [
            ROOT.kBlack,
            ROOT.kRed + 1,
            ROOT.kBlue + 1,
            ROOT.kGreen + 2,
            ROOT.kMagenta + 1,
            ROOT.kOrange + 7,
            ROOT.kCyan + 1,
        ]

    if line_styles is None:
        line_styles = [1] * len(hists)

    if marker_styles is None:
        marker_styles = [0] * len(hists)

    # Clone histograms so styling / normalization does not modify originals
    hists_draw = []
    for i, h in enumerate(hists):
        h_clone = h.Clone(f"{h.GetName()}_clone_{i}")
        h_clone.SetDirectory(0)

        if rebin:
            h_clone.Rebin(rebin)

        if normalize:
            integral = h_clone.Integral()
            if integral > 0:
                h_clone.Scale(1.0 / integral)

        if scale_factor:
            h_clone.Scale(scale_factor[i])

        color = colors[i % len(colors)]
        h_clone.SetLineColor(color)
        h_clone.SetMarkerColor(color)
        h_clone.SetLineWidth(line_width)
        h_clone.SetLineStyle(line_styles[i % len(line_styles)])
        h_clone.SetMarkerStyle(marker_styles[i % len(marker_styles)])
        h_clone.SetMarkerSize(marker_size)

        hists_draw.append(h_clone)

        if print_bins:
            for ib in range(0, h.GetNbinsX()+1):
                print(h.GetBinContent(ib))


    # Determine x range
    if x_range is not None:
        xmin, xmax = x_range
    else:
        xmin = min(h.GetXaxis().GetXmin() for h in hists_draw)
        xmax = max(h.GetXaxis().GetXmax() for h in hists_draw)

    # Determine y range
    if y_range is not None:
        ymin, ymax = y_range
    else:
        max_y = max(h.GetMaximum() for h in hists_draw)
        if logy:
            ymin = 0.5
            ymax = max_y * 50.0 if max_y > 0 else 10.0
        else:
            ymin = 0.0
            ymax = max_y * y_factor if max_y > 0 else 1.0

    # Create canvas
    cw, ch = canvas_size
    c = ROOT.TCanvas(f"c_{outname.replace('.', '_')}", "", cw, ch)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.16)
    c.SetRightMargin(0.05)

    if logy:
        c.SetLogy()

    if logx:
        c.SetLogx()

    # --- Handle bin labels if x_title is a list ---
    if isinstance(x_title, (list, tuple)):
        bin_labels = x_title

        n_bins = len(bin_labels)

        # Rebuild frame with correct binning
        h_frame = ROOT.TH1F(f"h_frame", "", n_bins, 0, n_bins)
        h_frame.SetDirectory(0)
        h_frame.SetMinimum(ymin)
        h_frame.SetMaximum(ymax)

        # Assign labels
        for i, lab in enumerate(bin_labels, start=1):
            h_frame.GetXaxis().SetBinLabel(i, lab)

        # Style for categorical axis
        h_frame.GetXaxis().LabelsOption("u")   # vertical labels (important!)
        
        h_frame.GetXaxis().SetLabelSize(0.045)
        h_frame.GetXaxis().SetTitle("")        # no numeric title
    else:
        # Standard numeric axis
        h_frame = ROOT.TH1F(f"h_frame", "", 1, xmin, xmax)
        h_frame.SetDirectory(0)
        h_frame.SetMinimum(ymin)
        h_frame.SetMaximum(ymax)

        if x_title == None:
            h_frame.GetXaxis().SetTitle(hists[0].GetXaxis().GetTitle())
        else:
            h_frame.GetXaxis().SetTitle(x_title)

    # Dummy histogram/frame
    # Use 1 bin only, just as axis frame
    ##h_frame = ROOT.TH1F("h_frame", "", 1, xmin, xmax)
    ##h_frame.SetMinimum(ymin)
    ##h_frame.SetMaximum(ymax)

    # Axis titles
    ##h_frame.GetXaxis().SetTitle(x_title)
    if y_title == None:
        h_frame.GetYaxis().SetTitle(hists[0].GetYaxis().GetTitle())
    else:
        h_frame.GetYaxis().SetTitle(y_title)

    # Make labels/titles large enough
    h_frame.GetXaxis().SetTitleSize(0.04)
    h_frame.GetYaxis().SetTitleSize(0.04)
    h_frame.GetXaxis().SetLabelSize(0.04)
    h_frame.GetYaxis().SetLabelSize(0.04)

    h_frame.GetXaxis().SetTitleOffset(1.10)
    h_frame.GetYaxis().SetTitleOffset(1.9)

    h_frame.GetXaxis().SetLabelOffset(0.01)
    h_frame.GetYaxis().SetLabelOffset(0.01)

    h_frame.GetXaxis().SetNdivisions(510)
    h_frame.GetYaxis().SetNdivisions(510)

    h_frame.Draw("AXIS")
    h_frame.Draw("AXIG SAME")  # helps make full frame look nice

    # Draw histograms
    first = True
    print(draw_option)
    for h in hists_draw:
        opt = draw_option
        if "text" in draw_option or "TEXT" in draw_option:
            ROOT.gStyle.SetPaintTextFormat(".2f")
            h.SetMarkerSize(1.5)
            h.Draw(f"{draw_option.replace('text', '').replace('TEXT', '')} SAME")
            h.Draw("TEXT SAME")
        else:
            h.Draw(f"{draw_option} SAME")
        first = False

    # Redraw axes on top
    h_frame.Draw("AXIS SAME")

    # Legend
    if len(labels) > 0:
        leg = make_dynamic_legend(
            hists_draw,
            labels,
            xmin,
            xmax,
            marker_size=marker_size,
            legend_header=legend_header,
            legend_ncols=legend_ncols
        )

        leg.Draw()

    # Optional text
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.16, 0.95, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.95, 0.95, extra_text_right)

    c.Modified()
    c.Update()
    c.SaveAs(f"{outname}.png")
    c.SaveAs(f"{outname}.pdf")


def plot_hist_2d(
    hist,
    outname,
    x_title="",
    y_title="",
    z_title="",
    x_range=None,
    y_range=None,
    z_range=None,
    logz=False,
    canvas_size=(900, 800),
    extra_text_left=None,
    extra_text_right=None,
    draw_option="COLZ",
):
    """
    Plot a single 2D histogram as a heat map.

    Parameters
    ----------
    hist : ROOT.TH2
        Input 2D histogram.
    outname : str
        Output file name, e.g. "plot.png".
    x_title, y_title, z_title : str
        Axis titles.
    x_range : tuple or None
        (xmin, xmax)
    y_range : tuple or None
        (ymin, ymax)
    z_range : tuple or None
        (zmin, zmax)
    logz : bool
        Use logarithmic z scale.
    canvas_size : tuple
        (width, height)
    extra_text : str or None
        Optional TLatex text.
    draw_option : str
        Usually "COLZ".
    """
    if not hist:
        raise ValueError("No histogram provided to plot_hist_2d().")

    h = hist.Clone(f"{hist.GetName()}_2d_clone")
    h.SetDirectory(0)

    # Apply displayed axis ranges if requested
    if x_range is not None:
        h.GetXaxis().SetRangeUser(x_range[0], x_range[1])
    if y_range is not None:
        h.GetYaxis().SetRangeUser(y_range[0], y_range[1])

    if z_range is not None:
        h.SetMinimum(z_range[0])
        h.SetMaximum(z_range[1])

    # Axis titles
    if x_title != None:
        h.GetXaxis().SetTitle(x_title)
    if y_title != None:
        h.GetYaxis().SetTitle(y_title)
    if z_title != None:
        h.GetZaxis().SetTitle(z_title)

    # Axis sizes
    h.GetXaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitleSize(0.04)
    h.GetZaxis().SetTitleSize(0.04)

    h.GetXaxis().SetLabelSize(0.04)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetZaxis().SetLabelSize(0.04)

    # Offsets
    h.GetXaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleOffset(1.35)
    h.GetZaxis().SetTitleOffset(1.20)

    h.GetXaxis().SetLabelOffset(0.01)
    h.GetYaxis().SetLabelOffset(0.01)
    h.GetZaxis().SetLabelOffset(0.01)

    h.GetXaxis().SetNdivisions(510)
    h.GetYaxis().SetNdivisions(510)
    h.GetZaxis().SetNdivisions(510)

    cw, ch = canvas_size
    c = ROOT.TCanvas(f"c_{h.GetName()}", "", cw, ch)

    # Slightly wider right margin for color palette
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.18)


    if logz:
        c.SetLogz()

    h.Draw(draw_option)

    # Update once so the palette exists
    c.Modified()
    c.Update()

    # Optional: style the palette axis a bit
    palette = h.GetListOfFunctions().FindObject("palette")
    if palette:
        palette.SetX1NDC(0.84)
        palette.SetX2NDC(0.89)

    # Optional text
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.14, 0.95, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.82, 0.95, extra_text_right)

    c.Modified()
    c.Update()
    c.SaveAs(f"{outname}.png")
    c.SaveAs(f"{outname}.pdf")