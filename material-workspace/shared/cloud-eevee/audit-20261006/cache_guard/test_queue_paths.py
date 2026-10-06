import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from queue_paths import execution_paths


class QueuePathTests(unittest.TestCase):
    def test_same_basename_in_two_runs_has_distinct_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            first = execution_paths(root / 'run_a/job.json', root / 'receipts')
            second = execution_paths(root / 'run_b/job.json', root / 'receipts')
            self.assertEqual(set(first), {'log', 'performance', 'cleanup'})
            for record in first:
                self.assertNotEqual(first[record], second[record])
                self.assertEqual(first[record].parent, root / 'receipts')
            self.assertEqual(len(set(first.values()) | set(second.values())), 6)

    def test_repeat_is_stable(self):
        self.assertEqual(execution_paths('/tmp/run/job.json', '/tmp/receipts'),
                         execution_paths('/tmp/run/job.json', '/tmp/receipts'))

    def test_equivalent_resolved_paths_are_the_same_job(self):
        self.assertEqual(execution_paths('/tmp/run/../run/job.json', '/tmp/receipts'),
                         execution_paths('/tmp/run/job.json', '/tmp/receipts'))

    def test_names_keep_readable_stem_and_short_digest(self):
        paths = execution_paths('/tmp/run/concrete.json', '/tmp/receipts')
        self.assertRegex(paths['log'].name, r'^concrete_[0-9a-f]{16}\.log$')
        self.assertRegex(paths['performance'].name, r'^concrete_[0-9a-f]{16}_performance\.json$')
        self.assertRegex(paths['cleanup'].name, r'^concrete_[0-9a-f]{16}_cleanup\.json$')


if __name__ == '__main__':
    unittest.main(verbosity=2)
