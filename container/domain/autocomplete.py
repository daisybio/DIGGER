from django.http import JsonResponse
from django.db import connection
from django.views.decorators.http import require_http_methods


# Derived from here https://stackoverflow.com/questions/32465052/using-typeahead-js-in-django-project
@require_http_methods(['GET'])
def gene_symbol_autocomplete(request):
    q = request.GET.get('q')
    o = request.GET.get('o')
    data = {}

    if q:
        # Define the raw SQL query
        raw_sql = """
        SELECT gene_symbol, ensembl_id
        FROM domain_gene_""" + o + """
        WHERE gene_symbol ILIKE %s
        LIMIT 20
        """

        try:
            # Execute the raw SQL query
            with connection.cursor() as cursor:
                cursor.execute(raw_sql, [q + '%'])  # Add the '%' wildcard to the query string
                rows = cursor.fetchall()

            # Format the results
            data = [{'gene_symbol': row[0], 'ensembl_id': row[1]} for row in rows]
        except IndexError:
            data = []

    return JsonResponse(data, safe=False)
