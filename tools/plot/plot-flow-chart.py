# -*- coding: utf-8 -*-
#@note: conda activate my-rdkit-env
#@note: https://www.rdkit.org/docs/GettingStartedInPython.html#drawing-molecules
#@note: https://www.rdkit.org/docs/Install.html#how-to-install-rdkit-with-conda

from graphviz import Digraph

#dot = Digraph(comment='workflow', format="png")
dot = Digraph(comment='workflow', format="svg")

dot.attr('node', shape='box')
dot.node('A', '整个纳米体系')

dot.node('B', '(嵌入)\n 分子片段一')
dot.node('C', '(嵌入)\n 分子片段二')
dot.node('D', '(嵌入)\n 分子片段三')

dot.node('E', '分子力学力场模型')
dot.node('F', '杂化模型')

dot.edges(['AB', 'AC', 'AD', 'BE', 'CE', 'DE', 'EF'])

# 在创建两圆点之间创建一条边
#dot.edge('B', 'C', 'test')


# 保存source到文件，并提供Graphviz引擎
dot.save('test-table.gv')  # 保存
dot.render('test-table.gv')
# dot.view()  # 显示

# 从保存的文件读取并显示

from graphviz import Source

s = Source.from_file('test-table.gv')
print(s.source)  # 打印代码
# s.view()  # 显示
