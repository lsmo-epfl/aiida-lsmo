#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt



'/Users/leopold/Personal/Postdoc-MARVEL/Projects/2017-09-15_pawel_plugin/zeoplusplus/network'  '-r' 'UFF.rad' '-ha' 'DEF' '-block' '1.525' '10' 'out.block' '-volpo' '1.525' '1.525' '1000' 'out.volpo' 'Mg_MOF_74.cif'   