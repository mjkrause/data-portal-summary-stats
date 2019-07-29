from chalice import Chalice
app = Chalice(app_name='portal-summary-stats-lambda')


@app.route('/')
def index():
    return {'sanity': 'OK'}

