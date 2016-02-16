import unittest,os,gzip

import splitter

class TestCase(unittest.TestCase):
   def setUp(self) :
      self.filename = os.path.join('.','test_data','chembl_export.txt.gz')
      self.lines = [x.decode('utf-8').strip().split() for x in gzip.open(self.filename,'r').readlines()]
      self.lines.pop(0)
      self.expr = splitter.coreExpr
   def test1Coverage(self) :
      " the pattern should consume the entire InChI "
      for line in self.lines:
          if len(line)<6:
              continue
          sinchi = line[1]
          inchi = line[3]
          smiles = line[5]
          #print(sinchi)
          match = self.expr.match(sinchi)
          self.assertTrue(match,smiles)
          self.assertEqual(match.span()[1],len(sinchi),'%s\n\n%s'%(smiles,match.string[match.end():]))
          #print(inchi)
          match = self.expr.match(inchi)
          self.assertTrue(match,smiles)
          self.assertEqual(match.span()[1],len(inchi),smiles)

   def _test2FullCoverage(self) :
      " the pattern should consume the entire InChI "
      filename = os.path.join('.','test_data','chembl_export.full.txt.gz')  # Note: this file isn't checked in because it's huge
      reader=gzip.open(filename,'r')
      for i,line in enumerate(reader):
          if not i:
              continue
          if not (i+1)%100000:
              print("Done:",i+1)
          line = line.decode('utf-8').strip().split()
          if len(line)<6:
              continue
          sinchi = line[1]
          inchi = line[3]
          smiles = line[5]

          match = self.expr.match(sinchi)
          self.assertTrue(match,'%s\n\n%s\n\n%s'%(smiles,match.string,match.string[match.end():]))
          self.assertEqual(match.span()[1],len(sinchi),'%s\n\n%s\n\n%s'%(smiles,match.string,match.string[match.end():]))
          match = self.expr.match(inchi)
          self.assertTrue(match,'%s\n\n%s\n\n%s'%(smiles,match.string,match.string[match.end():]))
          self.assertEqual(match.span()[1],len(inchi),'%s\n\n%s\n\n%s'%(smiles,match.string,match.string[match.end():]))


if __name__ == '__main__':
   unittest.main()
