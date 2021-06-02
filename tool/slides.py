import sys
import os

import numpy as np
import pandas as pd


def tikz_cells(_outfile, _data, _name, _prefix, _scale, _min_cell_width, _setting, _column_head_color,
               _extra_cell_content=r"\ldots"):
    _outfile.write(r"  %% Fill Color Styles" + "\n")
    _outfile.write(r"  \tikzstyle{CellHeaderFill} = [" + _column_head_color + "]" + "\n")
    _outfile.write(r"  \tikzstyle{CellFill} = [" + _setting["cell_color"] + "]" + "\n")
    _outfile.write(r"  \tikzstyle{CellHlFill} = [" + _setting["cell_color_hl"] + "]" + "\n")
    _outfile.write("\n")

    _outfile.write(r"  %% Element Styles" + "\n")
    _outfile.write(r"  \tikzstyle{CellHeader} = [CellHeaderFill, minimum width=" + _min_cell_width
                   + ", minimum height=1.5em, node distance=1.5em]" + "\n")
    _outfile.write(r"  \tikzstyle{Cell} = [CellFill, minimum width=" + _min_cell_width
                   + ", minimum height=1.5em, node distance=1.5em]" + "\n")
    _outfile.write(r"  \tikzstyle{CellHl} = [CellHlFill, minimum width=" + _min_cell_width
                   + ", minimum height=1.5em, node distance=1.5em]" + "\n")
    _outfile.write("\n")

    _outfile.write(r"  \begin{scope}[name prefix = " + _prefix + "]" + "\n")
    _outfile.write(r"    \node[name=name, CellHeader]{" + _name + "};" + "\n")
    _outfile.write(r"    \node[name=prev, below of=name, Cell] {" + _extra_cell_content + "};" + "\n")
    prev_cell = "prev"

    for i in range(0, len(_data)):
        cell_class = "CellHl" if (i in _setting["cell_hl"]) else "Cell"
        _outfile.write(r"    \node[name=" + str(i) + ", below of=" + prev_cell + ", " + cell_class + "] {"
                       + str(_data[i]) + "};" + "\n")
        prev_cell = str(i)

    _outfile.write(r"    \node[name=next, below of=" + prev_cell + r", Cell] {" + _extra_cell_content + "};" + "\n")

    _outfile.write(r"  \end{scope}" + "\n")


def tikz_column(_outfile, _data, _name, _prefix, _scale, _min_cell_width, _setting, _column_head_color):
    _outfile.write(r"\documentclass{standalone}" + "\n")
    _outfile.write(r"\usepackage{tikz}" + "\n")
    _outfile.write(r"\usetikzlibrary {positioning,shapes,calc,arrows.meta,graphs}" + "\n")
    _outfile.write("\n")
    _outfile.write(r"\begin{document}" + "\n")

    _outfile.write(r"\begin{tikzpicture}[font=\ttfamily \bfseries \small, scale=" + str(_scale)
                   + ", transform shape]" + "\n")

    tikz_cells(_outfile, _data, _name, _prefix + "-", _scale, _min_cell_width, _setting, _column_head_color)

    _outfile.write(r"\end{tikzpicture}" + "\n")
    _outfile.write(r"\end{document}" + "\n")


def create_tikz_column(_columns, _settings, _data, _scale):
    for l_column in _columns:
        print(l_column["name"])
        for l_setting in _settings:
            with open(output_dir + '/' + basename + '/' + l_column["file_name"] + l_setting["suffix_filename"] + ".tex",
                      'w') as l_outfile:
                tikz_column(l_outfile, _data[l_column["data_name"]], l_column["name"], l_column["data_name"], _scale,
                            str(1 / _scale) + r"\linewidth", l_setting, l_column["head_color"])


def tikz_cell_arcs(_outfile, _arcs, _prefix):
    _outfile.write(r"\tikzset{ >=stealth } %Define standard arrow tip" + "\n")

    for i in range(0, len(_arcs)):
        _outfile.write(r"\path[->] (" + _prefix + "-" + str(_arcs[i]["source"]) + ".west) edge[bend right] (" + _prefix
                       + "-" + (str(_arcs[i]["target"]) if _arcs[i]["target"] != -1 else "prev") + ".west);" + "\n")


def tikz_linked_column(_outfile, _data, _name, _prefix, _scale, _min_cell_width, _setting, _column_head_color, _arcs):
    _outfile.write(r"\documentclass{standalone}" + "\n")
    _outfile.write(r"\usepackage{tikz}" + "\n")
    _outfile.write(r"\usetikzlibrary {positioning,shapes,calc,arrows.meta,graphs}" + "\n")
    _outfile.write("\n")
    _outfile.write(r"\begin{document}" + "\n")

    _outfile.write(r"\begin{tikzpicture}[font=\ttfamily \bfseries \small, scale=" + str(_scale)
                   + ", transform shape]" + "\n")

    tikz_cells(_outfile, _data, _name, _prefix + "-", _scale, _min_cell_width, _setting, _column_head_color, "")
    _outfile.write("\n")
    tikz_cell_arcs(_outfile, _arcs, _prefix)

    _outfile.write(r"\end{tikzpicture}" + "\n")

    _outfile.write(r"\end{document}" + "\n")


def create_tikz_linked_column(_columns, _settings, _data, _scale):
    for l_column in _columns:
        print(l_column["name"])
        for l_setting in _settings:
            with open(output_dir + '/' + basename + '/' + l_column["file_name"] + l_setting["suffix_filename"] + ".tex",
                      'w') as l_outfile:
                # tikz_linked_column(l_outfile, _data[l_column["data_name"]], l_column["name"], l_column["data_name"],
                #                    _scale, str(1 / _scale * 0.7) + r"\linewidth", l_setting, l_column["head_color"],
                #                    l_setting["arcs"])
                tikz_linked_column(l_outfile, _data[l_column["data_name"]], l_column["name"], l_column["data_name"],
                                   _scale, "0", l_setting, l_column["head_color"],
                                   l_setting["arcs"])


def tikz_layer_cell(_outfile, _name, _min_layer_width, _layer_height, _layer_position, _layer_color):
    _outfile.write(r"  %% Fill Color Styles" + "\n")
    _outfile.write(r"  \tikzstyle{LayerFill} = [" + _layer_color + "]" + "\n")
    _outfile.write(r"  %% Element Styles" + "\n")
    # _outfile.write(r"  \tikzstyle{LayerHeader} = [minimum width=" + _min_layer_width
    #                + ", minimum height=1.5em, node distance=1.5em]" + "\n")
    _outfile.write(r"  \tikzstyle{Layer} = [LayerFill, minimum width=" + _min_layer_width
                   + ", minimum height=" + str(_layer_height) + "*1.5em, node distance=1.5em]" + "\n")

    _outfile.write("\n")
    _outfile.write(r"  \node[name=layer-body, " + _layer_position + ", Layer] {" + _name + "};" + "\n")


# def tikz_layer(_outfile, _name, _scale, _min_layer_width, _layer_height, _layer_color):
def tikz_layer(_outfile, _data, _name, _prefix, _scale, _min_cell_width, _setting, _column_head_color, _layer_name,
               _layer_height, _layer_color):
    _outfile.write(r"\documentclass{standalone}" + "\n")
    _outfile.write(r"\usepackage{tikz}" + "\n")
    _outfile.write(r"\usetikzlibrary {positioning,shapes,calc,arrows.meta,graphs}" + "\n")
    _outfile.write("\n")
    _outfile.write(r"\begin{document}" + "\n")

    _outfile.write(r"\begin{tikzpicture}[font=\ttfamily \bfseries \small, scale=" + str(_scale)
                   + ", transform shape]" + "\n")

    tikz_cells(_outfile, _data, _name, _prefix + "-", _scale, "0", _setting, _column_head_color, "")
    _outfile.write("\n")

    # position = "right=1cm of " + _prefix + "-" + str(1 + len(data) / 2 - _layer_height / 2)
    # position = "right=1cm of " + _prefix + "-" + str(int(len(data) / 2))
    cell_width = 0.7
    position = "right=" + str(1 - cell_width) + "*" + _min_cell_width + " of " + _prefix + "-" + str(int(len(data) / 2))
    tikz_layer_cell(_outfile, _layer_name, str(cell_width) + "*" + _min_cell_width, _layer_height, position,
                    _layer_color)
    _outfile.write("\n")

    _outfile.write(r"\path[-] (layer-body.north west) edge[dashed] (empty-prev.south west);" + "\n")
    _outfile.write(r"\path[-] (layer-body.south west) edge[dashed] (empty-next.north west);" + "\n")

    _outfile.write(r"\end{tikzpicture}" + "\n")

    _outfile.write(r"\end{document}" + "\n")


# def create_tikz_layer(_filename, _name, _scale, _n_data, _layer_color):
#     with open(output_dir + '/' + basename + '/' + _filename, 'w') as l_outfile:
#         tikz_layer(l_outfile, _name, _scale, str(1 / _scale) + r"\linewidth", _n_data + 3, _layer_color)


def create_tikz_layer(_columns, _settings, _data, _scale, _layer_color):
    for l_column in _columns:
        print(l_column["name"])
        for l_setting in _settings:
            with open(output_dir + '/' + basename + '/' + l_column["file_name"] + l_setting["suffix_filename"] + ".tex",
                      'w') as l_outfile:
                tikz_layer(l_outfile, _data[l_column["data_name"]], l_column["name"], l_column["data_name"], _scale,
                           str(1 / _scale) + r"\linewidth", l_setting, l_column["head_color"], l_column["layer_name"],
                           l_column["layer_height"], _layer_color)


def tikz_grammar_tree(_outfile, _data, _name, _prefix, _scale, _min_cell_width, _setting, _column_head_color, _grm,
                      _grm_sty):
    _outfile.write(r"\documentclass{standalone}" + "\n")
    _outfile.write(r"\usepackage{tikz}" + "\n")
    _outfile.write(r"\usetikzlibrary {positioning,shapes,calc,arrows.meta,graphs}" + "\n")
    _outfile.write("\n")
    _outfile.write(r"\begin{document}" + "\n")

    _outfile.write(r"\begin{tikzpicture}[font=\ttfamily \bfseries \small, scale=" + str(_scale)
                   + ", transform shape, new set=import nodes]" + "\n")

    tikz_cells(_outfile, _data, _name, _prefix + ", nodes={set=import nodes}", _scale, "1pt", _setting,
               _column_head_color, "")
    _outfile.write("\n")

    _outfile.write(_grm_sty)
    _outfile.write("\n")

    for l_node_id in _setting["node_hl"]:
        _outfile.write(r"\tikzset{" + l_node_id + "sty/.style={" + _setting["cell_color_hl"] + "}}" + "\n")

    # _outfile.write(r"\graph[grow left sep=1cm, "
    #                r"nodes={xshift=1\linewidth, yshift=-1cm, circle, minimum size=0.8cm,draw}, "
    #                r"edges={->,>=stealth}] {" + "\n")
    _outfile.write(r"\graph[grow left sep=.8cm, "
                   r"nodes={xshift=9cm, yshift=-1cm, circle, minimum size=0.8cm,draw,"
                   + _setting["node_opts"] + "}, "
                                             r"edges={->,>=stealth}] {" + "\n")

    _outfile.write(r"(import nodes);" + "\n")

    l_grm = _grm.replace("XXX", _prefix)
    _outfile.write(l_grm)

    _outfile.write(r"};" + "\n")

    _outfile.write(_setting["extra_cmp"] + "\n")

    _outfile.write(r"\end{tikzpicture}" + "\n")

    _outfile.write(r"\end{document}" + "\n")


def create_tikz_grammar_tree(_columns, _settings, _data, _scale, _grm, _grm_sty):
    for l_column in _columns:
        print(l_column["name"])
        for l_setting in _settings:
            with open(output_dir + '/' + basename + '/' + l_column["file_name"] + l_setting["suffix_filename"] + ".tex",
                      'w') as l_outfile:
                tikz_grammar_tree(l_outfile, _data[l_column["data_name"]], l_column["name"], l_column["data_name"],
                                  _scale, str(1 / _scale) + r"\linewidth", l_setting, l_column["head_color"], _grm,
                                  _grm_sty)


if __name__ == "__main__":
    # Count the arguments
    arguments = len(sys.argv) - 1

    if arguments < 3:
        print("The script must be called with: "
              "basename, "
              "doc, "
              "output directory")
        exit()

    basename = sys.argv[1]
    document = sys.argv[2]
    output_dir = sys.argv[3]

    # Data for collection
    data = pd.read_csv(basename + ".csv", delimiter=',', encoding='utf-8')
    data["suffix"] = [s + r"\ldots" for s in data["Suffix"]]
    data["empty"] = ["" for i in range(0, len(data.values))]

    info = pd.read_csv(basename + ".info.csv", delimiter=',', encoding='utf-8')

    if not os.path.exists(output_dir + '/' + basename):
        os.makedirs(output_dir + '/' + basename)

    column_scale = 0.45
    # cell_color = "fill=yellow!15"
    col_h_color = "fill=gray!15"
    cell_color_n = ""
    cell_color_hl = "fill=mLightGreen!35"
    cell_color_alert = "fill=mLightBrown!45"

    cell_range = {"begin": int(info["range_begin"]) - int(info["erange_begin"]),
                  "end": int(info["range_end"]) - int(info["erange_begin"])}
    range_s = range(cell_range["begin"], cell_range["end"] + 1)
    range_s_d = [i for i in range_s if data["DA"][i] == int(document)]

    start_pos = data["Position"][0]
    sada_arcs = [{"source": i,
                  "target": (data["SADA C Backward"][i] - start_pos) if data["SADA C Backward"][i] >= start_pos else -1}
                 for i in range(0, len(data.values))]
    sada_arcs_r = [a for a in sada_arcs if cell_range["begin"] <= a["source"] <= cell_range["end"]]
    sada_arcs_e = [a for a in sada_arcs if a["target"] < cell_range["begin"] <= a["source"]]

    range_s_e_l = [a["source"] for a in sada_arcs_e]

    range_bv = range(12, 19)

    # Basics
    settings = [
        {"suffix_filename": "", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_n, "cell_hl": []},
        {"suffix_filename": "-hl", "col_head_color": col_h_color, "cell_color": cell_color_hl,
         "cell_color_hl": cell_color_hl, "cell_hl": []},
        {"suffix_filename": "-r", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": range_s},
        {"suffix_filename": "-r-d", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": range_s_d},
        {"suffix_filename": "-r-d-e", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": [range_s_d[0], range_s_d[-1]]},
        {"suffix_filename": "-r-e-l", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": range_s_e_l},
        {"suffix_filename": "-bv", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": range_bv},
    ]

    columns = [
        {"file_name": "sa", "data_name": "SA", "name": "SA", "head_color": col_h_color},
        {"file_name": "da", "data_name": "DA", "name": "DA", "head_color": col_h_color},
        {"file_name": "suffix", "data_name": "suffix", "name": "Suffix", "head_color": col_h_color},
        # {"file_name": "csa", "data_name": "SA", "name": "CSA", "head_color": col_h_color},
        {"file_name": "csa", "data_name": "SA", "name": "$r$-index", "head_color": col_h_color},
        {"file_name": "pos", "data_name": "Position", "name": r"$\#$", "head_color": ""},
        {"file_name": "sada_c", "data_name": "SADA C Backward", "name": "C", "head_color": col_h_color},
        {"file_name": "ilcp", "data_name": "ILCP Backward", "name": "ILCP", "head_color": col_h_color},
    ]

    create_tikz_column(columns, settings, data, column_scale)

    # Data for individual documents
    data_doc = pd.read_csv(basename + "." + document + ".csv", delimiter=',', encoding='utf-8')
    range_d_s = range(3, len(data_doc.values) - 3)
    settings_doc = [
        {"suffix_filename": "", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": []},
        {"suffix_filename": "-r", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": range_d_s},
        {"suffix_filename": "-r-e", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": [range_d_s[0], range_d_s[-1]]},
        {"suffix_filename": "-r-e-alert", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_alert, "cell_hl": [range_d_s[0], range_d_s[-1]]},
    ]

    columns_doc = [
        # {"file_name": "csa_d", "data_name": "SA_d", "name": "CSA$_" + document + "$", "head_color": col_h_color},
        {"file_name": "csa_d", "data_name": "SA_d", "name": "$r$-index$_" + document + "$", "head_color": col_h_color},
        {"file_name": "pos_d", "data_name": "Position_d", "name": r"$\#$", "head_color": ""},
    ]

    create_tikz_column(columns_doc, settings_doc, data_doc, column_scale)

    # SADA
    settings_sada = [
        {"suffix_filename": "", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": [], "arcs": sada_arcs},
        # {"suffix_filename": "-hl", "col_head_color": col_h_color, "cell_color": cell_color_hl,
        #  "cell_color_hl": cell_color_hl, "cell_hl": [], "arcs": arcs},
        {"suffix_filename": "-r", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": [], "arcs": sada_arcs_r},
        {"suffix_filename": "-r-e", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": cell_color_hl, "cell_hl": [], "arcs": sada_arcs_e},
        # {"suffix_filename": "-r-e", "col_head_color": col_h_color, "cell_color": cell_color_n,
        #  "cell_color_hl": cell_color_n, "cell_hl": [], "arcs": []},
    ]

    columns_sada = [
        # {"file_name": "sada_c_arcs", "data_name": "SADA C Backward", "name": "C", "head_color": col_h_color},
        {"file_name": "sada_c_arcs", "data_name": "empty", "name": "", "head_color": ""},
    ]

    create_tikz_linked_column(columns_sada, settings_sada, data, column_scale)

    settings_sada_rmq = [settings[0]]

    columns_sada_rmq = [
        # {"file_name": "sada_c_arcs", "data_name": "SADA C Backward", "name": "C", "head_color": col_h_color},
        {"file_name": "sada_c_rmq", "data_name": "empty", "name": "", "head_color": "", "layer_name": r"RMQ$_{C}$",
         "layer_height": 15},
        {"file_name": "ilcp_rmq", "data_name": "empty", "name": "", "head_color": "", "layer_name": r"RMQ$_{ILCP}$",
         "layer_height": 15},
        {"file_name": "cilcp_rmq", "data_name": "empty", "name": "", "head_color": "",
         "layer_name": r"RMQ$_{ILCP^{\bigstar}}$", "layer_height": 10},
    ]
    # create_tikz_layer("sada_c_rmq.tex", r"RMQ$_{C}$", column_scale, len(data.values), col_h_color)
    create_tikz_layer(columns_sada_rmq, settings_sada_rmq, data, column_scale, col_h_color)

    # ILCP
    ilcp_run = [i for i in range(0, len(data["ILCP Run Heads"])) if data["ILCP Run Heads"][i] == 1]
    cilcp_run = [i for i in range(0, len(data["CILCP Run Heads"])) if data["CILCP Run Heads"][i] == 1]

    settings_ilcp_run = [
        {"suffix_filename": "ilcp_run", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": col_h_color + ",draw", "cell_hl": ilcp_run},
        {"suffix_filename": "cilcp_run", "col_head_color": col_h_color, "cell_color": cell_color_n,
         "cell_color_hl": col_h_color + ",draw", "cell_hl": cilcp_run},
    ]

    columns_ilcp_run = [
        {"file_name": "", "data_name": "empty", "name": "", "head_color": ""},
    ]

    create_tikz_column(columns_ilcp_run, settings_ilcp_run, data, column_scale)

    # GCDA
    gcda_data = pd.read_csv(basename + ".grm.csv", delimiter=',', encoding='utf-8', dtype=str)
    settings_gcda_grm = [
        {"suffix_filename": "gcda_grm", "col_head_color": col_h_color, "cell_color": cell_color_n + ", draw",
         "cell_color_hl": col_h_color + ",draw", "cell_hl": [], "node_opts": col_h_color, "node_hl": [],
         "extra_cmp": ""},
        {"suffix_filename": "gcda_grm-hl", "col_head_color": col_h_color, "cell_color": cell_color_n + ", draw",
         "cell_color_hl": col_h_color + ",draw", "cell_hl": [], "node_opts": cell_color_hl, "node_hl": [],
         "extra_cmp": ""},
        {"suffix_filename": "gcda_grm-bv", "col_head_color": col_h_color, "cell_color": cell_color_n + ", draw",
         "cell_color_hl": cell_color_hl + ",draw", "cell_hl": range_bv, "node_opts": col_h_color,
         "node_hl": ["18r17l14"],
         "extra_cmp": r"\node [circle split, draw, right=5cm of empty20, " + cell_color_hl
                      + r",align=center] { {\Large " + str(gcda_data["Var"][0]) + r"} \nodepart{lower} {\large "
                      + gcda_data["BV_F"][0] + r"} \\ {\large " + gcda_data["BV_L"][0] + "}};"},
        {"suffix_filename": "gcda_grm-r-c", "col_head_color": col_h_color, "cell_color": cell_color_n + ", draw",
         "cell_color_hl": cell_color_hl + ",draw", "cell_hl": [21], "node_opts": col_h_color,
         "node_hl": ["18l16l13r11r8", "18l16r12", "18r17r15l5l3", "18r17l14"],
         "extra_cmp": ""},
        {"suffix_filename": "gcda_grm-pdl", "col_head_color": col_h_color, "cell_color": cell_color_n + ", draw",
         "cell_color_hl": cell_color_hl + ",draw", "cell_hl": range_bv, "node_opts": col_h_color,
         "node_hl": ["18r17l14"],
         "extra_cmp": r"\node [circle split, draw, right=5cm of empty20, " + cell_color_hl
                      + r",align=center] { {\Large " + str(gcda_data["Var"][0])
                      + r"} \nodepart{lower} {\large 0:6} \\ {\large 2:1}};"},
    ]

    columns_gcda_grm = [
        {"file_name": "", "data_name": "empty", "name": "", "head_color": ""},
    ]

    with open(basename + ".grm", 'r') as grm_file:
        grm = grm_file.read()
    with open(basename + ".grm.sty", 'r') as grm_sty_file:
        grm_sty = grm_sty_file.read()
    create_tikz_grammar_tree(columns_gcda_grm, settings_gcda_grm, data, column_scale, grm, grm_sty)
