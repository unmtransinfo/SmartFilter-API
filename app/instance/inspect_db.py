import sqlite3, json, os

# Build absolute path to the correct database inside the instance folder
db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'pains.db'))
print(f"Opening DB at: {db_path}")

conn = sqlite3.connect(db_path)
cur = conn.cursor()

# Show tables to verify the table exists
cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cur.fetchall()
print("Tables in DB:", tables)

# Check if pains_run exists before querying it
if any("pains_run" in table for table in tables):
    cur.execute("SELECT id, timestamp, inputs, passed, failed FROM pains_run;")
    for id, ts, inputs_json, passed_json, failed_json in cur.fetchall():
        inputs = json.loads(inputs_json)
        passed = json.loads(passed_json)
        failed = json.loads(failed_json)
        print(f"Run {id} @ {ts}")
        print("  Inputs:", inputs)
        print("  Passed:", passed)
        print("  Failed:", ", ".join(failed))
        print()
else:
    print("Table 'pains_run' does not exist in this database.")

conn.close()
