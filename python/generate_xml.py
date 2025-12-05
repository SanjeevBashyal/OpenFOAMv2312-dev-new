import html

def create_node(id, value, x, y, width=160, height=60, style="rounded=1;whiteSpace=wrap;html=1;"):
    # Escape value for XML
    value = html.escape(value)
    return f'''    <mxCell id="{id}" value="{value}" style="{style}" vertex="1" parent="1">
      <mxGeometry x="{x}" y="{y}" width="{width}" height="{height}" as="geometry" />
    </mxCell>'''

def create_edge(id, source, target, text="", style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=0.5;exitY=1;exitDx=0;exitDy=0;entryX=0.5;entryY=0;entryDx=0;entryDy=0;"):
    return f'''    <mxCell id="{id}" value="{text}" style="{style}" edge="1" parent="1" source="{source}" target="{target}">
      <mxGeometry relative="1" as="geometry" />
    </mxCell>'''

def create_decision(id, value, x, y, width=100, height=100):
    style = "rhombus;whiteSpace=wrap;html=1;"
    return create_node(id, value, x, y, width, height, style)

nodes = []
edges = []

# Layout parameters
x_center = 400
y_start = 40
y_gap = 100

# 1. Start
nodes.append(create_node("start", "Start", x_center-40, y_start, 80, 40, "ellipse;whiteSpace=wrap;html=1;"))
y = y_start + 80

# 2. Initialization
nodes.append(create_node("init", "Initialize OpenFOAM\n(Mesh, Fields)\n&\nInitialize DEM\n(Particles, Walls, Contact)", x_center-100, y, 200, 80))
edges.append(create_edge("e1", "start", "init"))
y += 120

# 3. Settling Loop Start
nodes.append(create_decision("settling_cond", "Settling\nPhase?", x_center-50, y))
edges.append(create_edge("e2", "init", "settling_cond"))
y += 140

# 4. Settling Body
nodes.append(create_node("emit", "Emit Particles", x_center-80, y))
edges.append(create_edge("e3", "settling_cond", "emit", "Yes"))
y += 100

nodes.append(create_node("dem_forces", "Calc DEM Forces\n(Gravity + Contact)", x_center-80, y))
edges.append(create_edge("e4", "emit", "dem_forces"))
y += 100

nodes.append(create_node("dem_update", "Update Particles\n(Position, Velocity)", x_center-80, y))
edges.append(create_edge("e5", "dem_forces", "dem_update"))
y += 100

nodes.append(create_decision("settled", "Settled?", x_center-50, y))
edges.append(create_edge("e6", "dem_update", "settled"))

# Loop back for settling
# Edge from settled (No) to emit
# Need custom style for loop back
edges.append(f'''    <mxCell id="e_settle_loop" value="No" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=0;exitY=0.5;exitDx=0;exitDy=0;entryX=0;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="settled" target="emit">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
          <mxPoint x="{x_center-150}" y="{y+50}" />
          <mxPoint x="{x_center-150}" y="{y-200-100}" /> 
        </Array>
      </mxGeometry>
    </mxCell>''') # Approximate Y coords

y += 140

# 5. CFD Time Loop
nodes.append(create_decision("time_loop", "Run Time\nLoop?", x_center-50, y))
edges.append(create_edge("e7", "settled", "time_loop", "Yes"))
# Also edge from settling_cond (No) to time_loop
edges.append(f'''    <mxCell id="e_skip_settling" value="No" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=1;exitY=0.5;exitDx=0;exitDy=0;entryX=1;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="settling_cond" target="time_loop">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
          <mxPoint x="{x_center+150}" y="{y_start+200+50}" />
          <mxPoint x="{x_center+150}" y="{y+50}" />
        </Array>
      </mxGeometry>
    </mxCell>''')

y += 140

# 6. IBM
nodes.append(create_node("ibm", "Calculate IBM Fields\n(epsilon, Usolid)", x_center-80, y))
edges.append(create_edge("e8", "time_loop", "ibm", "Yes"))
y += 100

# 7. PIMPLE Loop
nodes.append(create_decision("pimple_loop", "PIMPLE\nLoop?", x_center-50, y))
edges.append(create_edge("e9", "ibm", "pimple_loop"))
y += 140

# 8. Fluid Solvers
nodes.append(create_node("alpha", "Solve Alpha Eqn\n(VOF / MULES)", x_center-80, y))
edges.append(create_edge("e10", "pimple_loop", "alpha", "Yes"))
y += 100

nodes.append(create_node("ueqn", "Solve Momentum Eqn\n(UEqn + IBM Penalty)", x_center-80, y))
edges.append(create_edge("e11", "alpha", "ueqn"))
y += 100

nodes.append(create_node("peqn", "Solve Pressure Eqn\n(pEqn + PIMPLE Corr)", x_center-80, y))
edges.append(create_edge("e12", "ueqn", "peqn"))
y += 100

nodes.append(create_node("turb", "Correct Turbulence", x_center-80, y))
edges.append(create_edge("e13", "peqn", "turb"))
y += 100

# 9. Coupling (Final Iter)
nodes.append(create_decision("final_iter", "Final\nIter?", x_center-50, y))
edges.append(create_edge("e14", "turb", "final_iter"))
y += 140

nodes.append(create_node("coupling", "Coupling:\n1. Calc Fluid Forces\n(Drag, Buoyancy)\n2. Add Gravity\n3. Update Particles", x_center-80, y, 200, 100))
edges.append(create_edge("e15", "final_iter", "coupling", "Yes"))

# Edge skipping coupling if not final
edges.append(f'''    <mxCell id="e_skip_coupling" value="No" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=0;exitY=0.5;exitDx=0;exitDy=0;entryX=0;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="final_iter" target="pimple_loop">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
           <mxPoint x="{x_center-200}" y="{y+50}" />
           <mxPoint x="{x_center-200}" y="{y-500-140}" /> 
        </Array>
      </mxGeometry>
    </mxCell>''')

y += 140

# Loop back PIMPLE
edges.append(f'''    <mxCell id="e_pimple_back" value="" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=0;exitY=0.5;exitDx=0;exitDy=0;entryX=0;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="coupling" target="pimple_loop">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
           <mxPoint x="{x_center-200}" y="{y+50}" />
           <mxPoint x="{x_center-200}" y="{y-640-140}" /> 
        </Array>
      </mxGeometry>
    </mxCell>''')

# 10. Write Data
nodes.append(create_node("write", "Write Data\n(Fields, Particles)", x_center-80, y))
# Edge from pimple_loop (No) to write
edges.append(f'''    <mxCell id="e_pimple_done" value="No" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=1;exitY=0.5;exitDx=0;exitDy=0;entryX=1;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="pimple_loop" target="write">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
          <mxPoint x="{x_center+200}" y="{y-640-140+50}" />
          <mxPoint x="{x_center+200}" y="{y+30}" />
        </Array>
      </mxGeometry>
    </mxCell>''')

y += 100

# Loop back Time
edges.append(f'''    <mxCell id="e_time_back" value="" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=0;exitY=0.5;exitDx=0;exitDy=0;entryX=0;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="write" target="time_loop">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
           <mxPoint x="{x_center-250}" y="{y-50+30}" />
           <mxPoint x="{x_center-250}" y="{y-1000-140}" /> 
        </Array>
      </mxGeometry>
    </mxCell>''')

# 11. End
nodes.append(create_node("end", "End", x_center-40, y, 80, 40, "ellipse;whiteSpace=wrap;html=1;"))
edges.append(f'''    <mxCell id="e_end" value="No" style="edgeStyle=orthogonalEdgeStyle;rounded=0;orthogonalLoop=1;jettySize=auto;html=1;exitX=1;exitY=0.5;exitDx=0;exitDy=0;entryX=1;entryY=0.5;entryDx=0;entryDy=0;" edge="1" parent="1" source="time_loop" target="end">
      <mxGeometry relative="1" as="geometry">
        <Array as="points">
          <mxPoint x="{x_center+250}" y="{y-1000-140+50}" />
          <mxPoint x="{x_center+250}" y="{y+20}" />
        </Array>
      </mxGeometry>
    </mxCell>''')


# Combine
xml_content = f'''<mxGraphModel dx="1213" dy="1500" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="850" pageHeight="2000" math="0" shadow="0">
  <root>
    <mxCell id="0" />
    <mxCell id="1" parent="0" />
{chr(10).join(nodes)}
{chr(10).join(edges)}
  </root>
</mxGraphModel>'''

print(xml_content)
