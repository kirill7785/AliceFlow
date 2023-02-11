object AddBlockForm: TAddBlockForm
  Left = 374
  Top = 167
  AutoSize = True
  Caption = 'Add Block'
  ClientHeight = 569
  ClientWidth = 369
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panelglobalconteiner: TPanel
    Left = 0
    Top = 0
    Width = 369
    Height = 569
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 0
    object RadioGroupglobalconteiner: TRadioGroup
      Left = 0
      Top = 13
      Width = 369
      Height = 52
      Caption = 'size and properties block'
      Columns = 4
      ItemIndex = 0
      Items.Strings = (
        'Info'
        'Geometry'
        'Properties'
        'Network')
      TabOrder = 0
      OnClick = RadioGroupglobalconteinerClick
    end
    object Panelinfo: TPanel
      Left = 9
      Top = 68
      Width = 352
      Height = 434
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 1
      object lblset_trans: TLabel
        Left = 11
        Top = 136
        Width = 84
        Height = 13
        Caption = 'Set Transparency'
      end
      object lbltransparent: TLabel
        Left = 16
        Top = 163
        Width = 57
        Height = 13
        Caption = 'Transparent'
      end
      object lblOpaque: TLabel
        Left = 88
        Top = 212
        Width = 38
        Height = 13
        Caption = 'Opaque'
      end
      object lblTransparentlab: TLabel
        Left = 176
        Top = 212
        Width = 57
        Height = 13
        Caption = 'Transparent'
      end
      object lbltransparencyvalue: TLabel
        Left = 160
        Top = 169
        Width = 6
        Height = 13
        Caption = '0'
      end
      object LabelLineWidth: TLabel
        Left = 24
        Top = 240
        Width = 48
        Height = 13
        Caption = 'LineWidth'
      end
      object GBname: TGroupBox
        Left = 8
        Top = 17
        Width = 313
        Height = 49
        Caption = 'block name'
        TabOrder = 0
        object Ename: TEdit
          Left = 8
          Top = 16
          Width = 289
          Height = 21
          TabOrder = 0
        end
      end
      object GroupBoxPriority: TGroupBox
        Left = 11
        Top = 71
        Width = 185
        Height = 49
        Caption = 'Priority'
        TabOrder = 1
        object EditPriority: TEdit
          Left = 13
          Top = 25
          Width = 156
          Height = 21
          TabOrder = 0
        end
      end
      object btncolor: TButton
        Left = 160
        Top = 126
        Width = 41
        Height = 25
        Caption = 'color'
        TabOrder = 2
        OnClick = btncolorClick
      end
      object scrlbrtrans1: TScrollBar
        Left = 89
        Top = 188
        Width = 145
        Height = 17
        PageSize = 0
        TabOrder = 3
        OnChange = scrlbrtrans1Change
      end
      object CheckBoxVisible: TCheckBox
        Left = 224
        Top = 96
        Width = 97
        Height = 17
        Caption = 'Visible'
        Checked = True
        State = cbChecked
        TabOrder = 4
      end
      object ComboBoxLineWidth: TComboBox
        Left = 88
        Top = 237
        Width = 36
        Height = 21
        ItemIndex = 0
        TabOrder = 5
        Text = '1'
        Items.Strings = (
          '1'
          '2'
          '3'
          '4'
          '5'
          '6')
      end
      object GroupBoxActivityCondition: TGroupBox
        Left = 16
        Top = 279
        Width = 305
        Height = 109
        Caption = 'Activity condition'
        TabOrder = 6
        object LabelOperation: TLabel
          Left = 128
          Top = 24
          Width = 46
          Height = 13
          Caption = 'Operation'
        end
        object LabelLeftExpression: TLabel
          Left = 16
          Top = 24
          Width = 74
          Height = 13
          Caption = ' Left expression'
        end
        object LabelRightExpression: TLabel
          Left = 204
          Top = 24
          Width = 78
          Height = 13
          Caption = 'Right expression'
        end
        object ComboBoxZnak: TComboBox
          Left = 128
          Top = 56
          Width = 49
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '=='
          Items.Strings = (
            '=='
            '!='
            '<'
            '<='
            '>'
            '>=')
        end
        object ComboBoxLeftExpression: TComboBox
          Left = 16
          Top = 43
          Width = 106
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Expression'
          OnChange = ComboBoxLeftExpressionChange
          Items.Strings = (
            'Expression'
            'obj X Min'
            'obj X Max'
            'obj Y Min'
            'obj Y Max'
            'obj Z Min'
            'obj Z Max')
        end
        object EditLeftExpression: TEdit
          Left = 16
          Top = 70
          Width = 105
          Height = 21
          TabOrder = 2
          Text = '1'
        end
        object EditRightExpression: TEdit
          Left = 183
          Top = 56
          Width = 106
          Height = 21
          TabOrder = 3
          Text = '1'
        end
      end
    end
    object Bapply: TButton
      Left = 272
      Top = 522
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = BapplyClick
    end
    object PanelGeometry: TPanel
      Left = 9
      Top = 71
      Width = 352
      Height = 434
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 3
      Visible = False
      object LabelGeometryType: TLabel
        Left = 32
        Top = 24
        Width = 66
        Height = 13
        Caption = 'geometry type'
      end
      object LabelPlane: TLabel
        Left = 30
        Top = 52
        Width = 27
        Height = 13
        Caption = 'Plane'
        Visible = False
      end
      object GBsizeBlock: TGroupBox
        Left = 20
        Top = 84
        Width = 312
        Height = 141
        Caption = 'size block'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        object LxS: TLabel
          Left = 8
          Top = 16
          Width = 12
          Height = 13
          Caption = 'xS'
        end
        object LyS: TLabel
          Left = 8
          Top = 48
          Width = 12
          Height = 13
          Caption = 'yS'
        end
        object LzS: TLabel
          Left = 8
          Top = 80
          Width = 12
          Height = 13
          Caption = 'zS'
        end
        object LxE: TLabel
          Left = 134
          Top = 23
          Width = 12
          Height = 13
          Caption = 'xE'
        end
        object LyE: TLabel
          Left = 134
          Top = 48
          Width = 12
          Height = 13
          Caption = 'yE'
        end
        object LzE: TLabel
          Left = 134
          Top = 85
          Width = 12
          Height = 13
          Caption = 'zE'
        end
        object Label1: TLabel
          Left = 104
          Top = 16
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label2: TLabel
          Left = 104
          Top = 48
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label3: TLabel
          Left = 104
          Top = 80
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label4: TLabel
          Left = 276
          Top = 24
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label5: TLabel
          Left = 274
          Top = 61
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label6: TLabel
          Left = 274
          Top = 80
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object ExS: TEdit
          Left = 32
          Top = 16
          Width = 57
          Height = 21
          TabOrder = 0
        end
        object EyS: TEdit
          Left = 32
          Top = 48
          Width = 57
          Height = 21
          TabOrder = 1
        end
        object EzS: TEdit
          Left = 32
          Top = 77
          Width = 57
          Height = 21
          TabOrder = 2
        end
        object ExE: TEdit
          Left = 205
          Top = 16
          Width = 65
          Height = 21
          TabOrder = 3
        end
        object EyE: TEdit
          Left = 205
          Top = 43
          Width = 65
          Height = 21
          TabOrder = 4
        end
        object EzE: TEdit
          Left = 205
          Top = 80
          Width = 65
          Height = 21
          TabOrder = 5
        end
        object CheckBoxCylinder2Prism: TCheckBox
          Left = 24
          Top = 112
          Width = 96
          Height = 17
          Caption = 'Cylinder2Prism'
          TabOrder = 6
        end
        object CheckBoxFixedCylinder: TCheckBox
          Left = 140
          Top = 107
          Width = 97
          Height = 17
          Caption = 'Fixed Cylinder'
          TabOrder = 7
        end
      end
      object ComboBoxgeometrytype: TComboBox
        Left = 112
        Top = 16
        Width = 193
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = 'Prism'
        OnChange = ComboBoxgeometrytypeChange
        Items.Strings = (
          'Prism'
          'Cylinder'
          'Polygon'
          'CAD')
      end
      object ComboBoxPlane: TComboBox
        Left = 72
        Top = 44
        Width = 54
        Height = 21
        ItemIndex = 0
        TabOrder = 2
        Text = 'X-Y'
        Visible = False
        Items.Strings = (
          'X-Y'
          'X-Z'
          'Y-Z')
      end
      object GBPolygonGeom: TGroupBox
        Left = 18
        Top = 84
        Width = 314
        Height = 157
        Margins.Left = 2
        Margins.Top = 2
        Margins.Right = 2
        Margins.Bottom = 2
        Caption = 'Location'
        TabOrder = 3
        Visible = False
        object Labeldim1: TLabel
          Left = 139
          Top = 27
          Width = 48
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labeldim1'
        end
        object Labeldimx: TLabel
          Left = 139
          Top = 53
          Width = 47
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labeldimx'
        end
        object Labeldimy: TLabel
          Left = 139
          Top = 75
          Width = 47
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labeldimy'
        end
        object Labeldimz: TLabel
          Left = 139
          Top = 99
          Width = 47
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labeldimz'
        end
        object LabelnameH: TLabel
          Left = 11
          Top = 32
          Width = 31
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Height'
        end
        object Labelx: TLabel
          Left = 11
          Top = 53
          Width = 31
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labelx'
        end
        object Labely: TLabel
          Left = 11
          Top = 80
          Width = 31
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labely'
        end
        object Labelz: TLabel
          Left = 13
          Top = 107
          Width = 31
          Height = 13
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Labelz'
        end
        object ListBoxvert: TListBox
          Left = 168
          Top = 16
          Width = 137
          Height = 107
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          ItemHeight = 13
          TabOrder = 0
          OnClick = ListBoxvertClick
        end
        object ButtonAdd: TButton
          Left = 176
          Top = 133
          Width = 50
          Height = 17
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Add'
          TabOrder = 1
          OnClick = ButtonAddClick
        end
        object ButtonRemove: TButton
          Left = 240
          Top = 133
          Width = 50
          Height = 17
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          Caption = 'Remove'
          TabOrder = 2
          OnClick = ButtonRemoveClick
        end
        object EditHeight: TEdit
          Left = 48
          Top = 27
          Width = 81
          Height = 21
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          TabOrder = 3
          Text = 'EditHeight'
        end
        object Editx: TEdit
          Left = 50
          Top = 51
          Width = 81
          Height = 21
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          TabOrder = 4
          Text = 'Editx'
        end
        object Edity: TEdit
          Left = 50
          Top = 75
          Width = 81
          Height = 21
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          TabOrder = 5
          Text = 'Edity'
        end
        object Editz: TEdit
          Left = 48
          Top = 99
          Width = 81
          Height = 21
          Margins.Left = 2
          Margins.Top = 2
          Margins.Right = 2
          Margins.Bottom = 2
          TabOrder = 6
          Text = 'Editz'
        end
      end
    end
    object PanelProperties: TPanel
      Left = 9
      Top = 71
      Width = 352
      Height = 434
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 4
      Visible = False
      object RadioGroupType: TRadioGroup
        Left = 10
        Top = 16
        Width = 289
        Height = 51
        Caption = 'Block type'
        Color = clMoneyGreen
        Columns = 3
        ItemIndex = 0
        Items.Strings = (
          'SOLID'
          'HOLLOW'
          'FLUID')
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        OnClick = RadioGroupTypeClick
      end
      object Panel1: TPanel
        Left = 9
        Top = 79
        Width = 288
        Height = 78
        Color = clMoneyGreen
        ParentBackground = False
        TabOrder = 1
        object Label7: TLabel
          Left = 16
          Top = 16
          Width = 99
          Height = 13
          Caption = 'Surface specification'
        end
        object ButtonRadiation: TButton
          Left = 13
          Top = 35
          Width = 75
          Height = 25
          Caption = 'Radiation'
          TabOrder = 0
          OnClick = ButtonRadiationClick
        end
      end
      object GroupBoxPropBl: TGroupBox
        Left = 11
        Top = 169
        Width = 313
        Height = 184
        Caption = 'Thermal specification'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 2
        object LMN: TLabel
          Left = 16
          Top = 48
          Width = 3
          Height = 13
          Color = clMoneyGreen
          ParentColor = False
        end
        object lblmateril: TLabel
          Left = 16
          Top = 29
          Width = 71
          Height = 13
          Caption = 'material name :'
        end
        object Labelpowerinfo: TLabel
          Left = 143
          Top = 42
          Width = 74
          Height = 13
          Caption = 'LabelPowerInfo'
        end
        object GBmaterial: TGroupBox
          Left = 13
          Top = 76
          Width = 289
          Height = 105
          Caption = 'Material'
          Color = clMoneyGreen
          ParentBackground = False
          ParentColor = False
          TabOrder = 0
          object RGSelect: TRadioGroup
            Left = 7
            Top = 16
            Width = 226
            Height = 49
            Caption = 'Define '
            Color = clMoneyGreen
            Columns = 2
            ItemIndex = 0
            Items.Strings = (
              'Program Library'
              'User-Defined')
            ParentColor = False
            TabOrder = 0
            OnClick = RGSelectClick
          end
          object BEditApply: TButton
            Left = 239
            Top = 21
            Width = 41
            Height = 81
            Caption = 'Edit'
            TabOrder = 1
            OnClick = BEditApplyClick
          end
          object CBselectAction: TComboBox
            Left = 8
            Top = 72
            Width = 225
            Height = 21
            TabOrder = 2
            Text = 'Edit Current Material'
            Items.Strings = (
              'Edit Current Material'
              'Create New Material'
              'Select Project Material ')
          end
        end
        object GBPower: TGroupBox
          Left = 231
          Top = 21
          Width = 65
          Height = 49
          Caption = 'Power'
          Color = clMoneyGreen
          ParentBackground = False
          ParentColor = False
          TabOrder = 1
          object ButtonTransient: TButton
            Left = -2
            Top = 21
            Width = 54
            Height = 25
            Caption = 'Edit'
            TabOrder = 0
            OnClick = ButtonTransientClick
          end
        end
      end
    end
    object PanelNetwork: TPanel
      Left = 8
      Top = 71
      Width = 353
      Height = 433
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 5
      Visible = False
      object GroupBoxinx: TGroupBox
        Left = 16
        Top = 27
        Width = 274
        Height = 62
        Caption = 'inx'
        TabOrder = 0
        object ComboBoxinx: TComboBox
          Left = 19
          Top = 25
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
            '1'
            '2'
            '3'
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33'
            '34'
            '35'
            '36'
            '37'
            '38'
            '39'
            '40'
            '41'
            '42'
            '43'
            '44'
            '45'
            '46'
            '47'
            '48'
            '49'
            '50'
            '51'
            '52'
            '53'
            '54'
            '55'
            '56'
            '57'
            '58'
            '59'
            '60'
            '61'
            '62'
            '63'
            '64'
            '65'
            '66'
            '67'
            '68'
            '69'
            '70'
            '71'
            '72'
            '73'
            '74'
            '75'
            '76'
            '77'
            '78'
            '79'
            '80'
            '81'
            '82'
            '83'
            '84'
            '85'
            '86'
            '87'
            '88'
            '89'
            '90'
            '91'
            '92'
            '93'
            '94'
            '95'
            '96'
            '97'
            '98'
            '99'
            '100'
            '101'
            '102'
            '103'
            '104'
            '105'
            '106'
            '107'
            '108'
            '109'
            '110'
            '111'
            '112'
            '113'
            '114'
            '115'
            '116'
            '117'
            '118'
            '119'
            '120'
            '121'
            '122'
            '123'
            '124'
            '125'
            '126'
            '127'
            '128'
            '129'
            '130'
            '131'
            '132'
            '133'
            '134'
            '135'
            '136'
            '137'
            '138'
            '139'
            '140'
            '141'
            '142'
            '143'
            '144'
            '145'
            '146'
            '147'
            '148'
            '149'
            '150'
            '151'
            '152'
            '153'
            '154'
            '155'
            '156'
            '157'
            '158'
            '159'
            '160'
            '161'
            '162'
            '163'
            '164'
            '165'
            '166'
            '167'
            '168'
            '169'
            '170'
            '171'
            '172'
            '173'
            '174'
            '175'
            '176'
            '177'
            '178'
            '179'
            '180'
            '181'
            '182'
            '183'
            '184'
            '185'
            '186'
            '187'
            '188'
            '189'
            '190'
            '191'
            '192'
            '193'
            '194'
            '195'
            '196'
            '197'
            '198'
            '199'
            '200'
            '201'
            '202'
            '203'
            '204'
            '205'
            '206'
            '207'
            '208'
            '209'
            '210'
            '211'
            '212'
            '213'
            '214'
            '215'
            '216'
            '217'
            '218'
            '219'
            '220'
            '221'
            '222'
            '223'
            '224'
            '225'
            '226'
            '227'
            '228'
            '229'
            '230'
            '231'
            '232'
            '233'
            '234'
            '235'
            '236'
            '237'
            '238'
            '239'
            '240'
            '241'
            '242'
            '243'
            '244'
            '245'
            '246'
            '247'
            '248'
            '249'
            '250'
            '251'
            '252'
            '253'
            '254'
            '255'
            '256'
            '257'
            '258'
            '259'
            '260'
            '261'
            '262'
            '263'
            '264'
            '265'
            '266'
            '267'
            '268'
            '269'
            '270'
            '271'
            '272'
            '273'
            '274'
            '275'
            '276'
            '277'
            '278'
            '279'
            '280'
            '281'
            '282'
            '283'
            '284'
            '285'
            '286'
            '287'
            '288'
            '289'
            '290'
            '291'
            '292'
            '293'
            '294'
            '295'
            '296'
            '297'
            '298'
            '299'
            '300'
            '301'
            '302'
            '303'
            '304'
            '305'
            '306'
            '307'
            '308'
            '309'
            '310'
            '311'
            '312'
            '313'
            '314'
            '315'
            '316'
            '317'
            '318'
            '319'
            '320'
            '321'
            '322'
            '323'
            '324'
            '325'
            '326'
            '327'
            '328'
            '329'
            '330'
            '331'
            '332'
            '333'
            '334'
            '335'
            '336'
            '337'
            '338'
            '339'
            '340'
            '341'
            '342'
            '343'
            '344'
            '345'
            '346'
            '347'
            '348'
            '349'
            '350'
            '351'
            '352'
            '353'
            '354'
            '355'
            '356'
            '357'
            '358'
            '359'
            '360'
            '361'
            '362'
            '363'
            '364'
            '365'
            '366'
            '367'
            '368'
            '369'
            '370'
            '371'
            '372'
            '373'
            '374'
            '375'
            '376'
            '377'
            '378'
            '379'
            '380'
            '381'
            '382'
            '383'
            '384'
            '385'
            '386'
            '387'
            '388'
            '389'
            '390'
            '391'
            '392'
            '393'
            '394'
            '395'
            '396'
            '397'
            '398'
            '399'
            '400'
            '401'
            '402'
            '403'
            '404'
            '405'
            '406'
            '407'
            '408'
            '409'
            '410'
            '411'
            '412'
            '413'
            '414'
            '415'
            '416'
            '417'
            '418'
            '419'
            '420'
            '421'
            '422'
            '423'
            '424'
            '425'
            '426'
            '427'
            '428'
            '429'
            '430'
            '431'
            '432'
            '433'
            '434'
            '435'
            '436'
            '437'
            '438'
            '439'
            '440'
            '441'
            '442'
            '443'
            '444'
            '445'
            '446'
            '447'
            '448'
            '449'
            '450'
            '451'
            '452'
            '453'
            '454'
            '455'
            '456'
            '457'
            '458'
            '459'
            '460'
            '461'
            '462'
            '463'
            '464'
            '465'
            '466'
            '467'
            '468'
            '469'
            '470'
            '471'
            '472'
            '473'
            '474'
            '475'
            '476'
            '477'
            '478'
            '479'
            '480'
            '481'
            '482'
            '483'
            '484'
            '485'
            '486'
            '487'
            '488'
            '489'
            '490'
            '491'
            '492'
            '493'
            '494'
            '495'
            '496'
            '497'
            '498'
            '499'
            '500')
        end
      end
      object GroupBoxiny: TGroupBox
        Left = 14
        Top = 96
        Width = 274
        Height = 67
        Caption = 'iny'
        TabOrder = 1
        object ComboBoxiny: TComboBox
          Left = 22
          Top = 27
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
            '1'
            '2'
            '3'
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33'
            '34'
            '35'
            '36'
            '37'
            '38'
            '39'
            '40'
            '41'
            '42'
            '43'
            '44'
            '45'
            '46'
            '47'
            '48'
            '49'
            '50'
            '51'
            '52'
            '53'
            '54'
            '55'
            '56'
            '57'
            '58'
            '59'
            '60'
            '61'
            '62'
            '63'
            '64'
            '65'
            '66'
            '67'
            '68'
            '69'
            '70'
            '71'
            '72'
            '73'
            '74'
            '75'
            '76'
            '77'
            '78'
            '79'
            '80'
            '81'
            '82'
            '83'
            '84'
            '85'
            '86'
            '87'
            '88'
            '89'
            '90'
            '91'
            '92'
            '93'
            '94'
            '95'
            '96'
            '97'
            '98'
            '99'
            '100'
            '101'
            '102'
            '103'
            '104'
            '105'
            '106'
            '107'
            '108'
            '109'
            '110'
            '111'
            '112'
            '113'
            '114'
            '115'
            '116'
            '117'
            '118'
            '119'
            '120'
            '121'
            '122'
            '123'
            '124'
            '125'
            '126'
            '127'
            '128'
            '129'
            '130'
            '131'
            '132'
            '133'
            '134'
            '135'
            '136'
            '137'
            '138'
            '139'
            '140'
            '141'
            '142'
            '143'
            '144'
            '145'
            '146'
            '147'
            '148'
            '149'
            '150'
            '151'
            '152'
            '153'
            '154'
            '155'
            '156'
            '157'
            '158'
            '159'
            '160'
            '161'
            '162'
            '163'
            '164'
            '165'
            '166'
            '167'
            '168'
            '169'
            '170'
            '171'
            '172'
            '173'
            '174'
            '175'
            '176'
            '177'
            '178'
            '179'
            '180'
            '181'
            '182'
            '183'
            '184'
            '185'
            '186'
            '187'
            '188'
            '189'
            '190'
            '191'
            '192'
            '193'
            '194'
            '195'
            '196'
            '197'
            '198'
            '199'
            '200'
            '201'
            '202'
            '203'
            '204'
            '205'
            '206'
            '207'
            '208'
            '209'
            '210'
            '211'
            '212'
            '213'
            '214'
            '215'
            '216'
            '217'
            '218'
            '219'
            '220'
            '221'
            '222'
            '223'
            '224'
            '225'
            '226'
            '227'
            '228'
            '229'
            '230'
            '231'
            '232'
            '233'
            '234'
            '235'
            '236'
            '237'
            '238'
            '239'
            '240'
            '241'
            '242'
            '243'
            '244'
            '245'
            '246'
            '247'
            '248'
            '249'
            '250'
            '251'
            '252'
            '253'
            '254'
            '255'
            '256'
            '257'
            '258'
            '259'
            '260'
            '261'
            '262'
            '263'
            '264'
            '265'
            '266'
            '267'
            '268'
            '269'
            '270'
            '271'
            '272'
            '273'
            '274'
            '275'
            '276'
            '277'
            '278'
            '279'
            '280'
            '281'
            '282'
            '283'
            '284'
            '285'
            '286'
            '287'
            '288'
            '289'
            '290'
            '291'
            '292'
            '293'
            '294'
            '295'
            '296'
            '297'
            '298'
            '299'
            '300'
            '301'
            '302'
            '303'
            '304'
            '305'
            '306'
            '307'
            '308'
            '309'
            '310'
            '311'
            '312'
            '313'
            '314'
            '315'
            '316'
            '317'
            '318'
            '319'
            '320'
            '321'
            '322'
            '323'
            '324'
            '325'
            '326'
            '327'
            '328'
            '329'
            '330'
            '331'
            '332'
            '333'
            '334'
            '335'
            '336'
            '337'
            '338'
            '339'
            '340'
            '341'
            '342'
            '343'
            '344'
            '345'
            '346'
            '347'
            '348'
            '349'
            '350'
            '351'
            '352'
            '353'
            '354'
            '355'
            '356'
            '357'
            '358'
            '359'
            '360'
            '361'
            '362'
            '363'
            '364'
            '365'
            '366'
            '367'
            '368'
            '369'
            '370'
            '371'
            '372'
            '373'
            '374'
            '375'
            '376'
            '377'
            '378'
            '379'
            '380'
            '381'
            '382'
            '383'
            '384'
            '385'
            '386'
            '387'
            '388'
            '389'
            '390'
            '391'
            '392'
            '393'
            '394'
            '395'
            '396'
            '397'
            '398'
            '399'
            '400'
            '401'
            '402'
            '403'
            '404'
            '405'
            '406'
            '407'
            '408'
            '409'
            '410'
            '411'
            '412'
            '413'
            '414'
            '415'
            '416'
            '417'
            '418'
            '419'
            '420'
            '421'
            '422'
            '423'
            '424'
            '425'
            '426'
            '427'
            '428'
            '429'
            '430'
            '431'
            '432'
            '433'
            '434'
            '435'
            '436'
            '437'
            '438'
            '439'
            '440'
            '441'
            '442'
            '443'
            '444'
            '445'
            '446'
            '447'
            '448'
            '449'
            '450'
            '451'
            '452'
            '453'
            '454'
            '455'
            '456'
            '457'
            '458'
            '459'
            '460'
            '461'
            '462'
            '463'
            '464'
            '465'
            '466'
            '467'
            '468'
            '469'
            '470'
            '471'
            '472'
            '473'
            '474'
            '475'
            '476'
            '477'
            '478'
            '479'
            '480'
            '481'
            '482'
            '483'
            '484'
            '485'
            '486'
            '487'
            '488'
            '489'
            '490'
            '491'
            '492'
            '493'
            '494'
            '495'
            '496'
            '497'
            '498'
            '499'
            '500')
        end
      end
      object GroupBoxinz: TGroupBox
        Left = 14
        Top = 169
        Width = 274
        Height = 70
        Caption = 'inz'
        TabOrder = 2
        object ComboBoxinz: TComboBox
          Left = 23
          Top = 32
          Width = 145
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = '1'
          Items.Strings = (
            '1'
            '2'
            '3'
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33'
            '34'
            '35'
            '36'
            '37'
            '38'
            '39'
            '40'
            '41'
            '42'
            '43'
            '44'
            '45'
            '46'
            '47'
            '48'
            '49'
            '50'
            '51'
            '52'
            '53'
            '54'
            '55'
            '56'
            '57'
            '58'
            '59'
            '60'
            '61'
            '62'
            '63'
            '64'
            '65'
            '66'
            '67'
            '68'
            '69'
            '70'
            '71'
            '72'
            '73'
            '74'
            '75'
            '76'
            '77'
            '78'
            '79'
            '80'
            '81'
            '82'
            '83'
            '84'
            '85'
            '86'
            '87'
            '88'
            '89'
            '90'
            '91'
            '92'
            '93'
            '94'
            '95'
            '96'
            '97'
            '98'
            '99'
            '100'
            '101'
            '102'
            '103'
            '104'
            '105'
            '106'
            '107'
            '108'
            '109'
            '110'
            '111'
            '112'
            '113'
            '114'
            '115'
            '116'
            '117'
            '118'
            '119'
            '120'
            '121'
            '122'
            '123'
            '124'
            '125'
            '126'
            '127'
            '128'
            '129'
            '130'
            '131'
            '132'
            '133'
            '134'
            '135'
            '136'
            '137'
            '138'
            '139'
            '140'
            '141'
            '142'
            '143'
            '144'
            '145'
            '146'
            '147'
            '148'
            '149'
            '150'
            '151'
            '152'
            '153'
            '154'
            '155'
            '156'
            '157'
            '158'
            '159'
            '160'
            '161'
            '162'
            '163'
            '164'
            '165'
            '166'
            '167'
            '168'
            '169'
            '170'
            '171'
            '172'
            '173'
            '174'
            '175'
            '176'
            '177'
            '178'
            '179'
            '180'
            '181'
            '182'
            '183'
            '184'
            '185'
            '186'
            '187'
            '188'
            '189'
            '190'
            '191'
            '192'
            '193'
            '194'
            '195'
            '196'
            '197'
            '198'
            '199'
            '200'
            '201'
            '202'
            '203'
            '204'
            '205'
            '206'
            '207'
            '208'
            '209'
            '210'
            '211'
            '212'
            '213'
            '214'
            '215'
            '216'
            '217'
            '218'
            '219'
            '220'
            '221'
            '222'
            '223'
            '224'
            '225'
            '226'
            '227'
            '228'
            '229'
            '230'
            '231'
            '232'
            '233'
            '234'
            '235'
            '236'
            '237'
            '238'
            '239'
            '240'
            '241'
            '242'
            '243'
            '244'
            '245'
            '246'
            '247'
            '248'
            '249'
            '250'
            '251'
            '252'
            '253'
            '254'
            '255'
            '256'
            '257'
            '258'
            '259'
            '260'
            '261'
            '262'
            '263'
            '264'
            '265'
            '266'
            '267'
            '268'
            '269'
            '270'
            '271'
            '272'
            '273'
            '274'
            '275'
            '276'
            '277'
            '278'
            '279'
            '280'
            '281'
            '282'
            '283'
            '284'
            '285'
            '286'
            '287'
            '288'
            '289'
            '290'
            '291'
            '292'
            '293'
            '294'
            '295'
            '296'
            '297'
            '298'
            '299'
            '300'
            '301'
            '302'
            '303'
            '304'
            '305'
            '306'
            '307'
            '308'
            '309'
            '310'
            '311'
            '312'
            '313'
            '314'
            '315'
            '316'
            '317'
            '318'
            '319'
            '320'
            '321'
            '322'
            '323'
            '324'
            '325'
            '326'
            '327'
            '328'
            '329'
            '330'
            '331'
            '332'
            '333'
            '334'
            '335'
            '336'
            '337'
            '338'
            '339'
            '340'
            '341'
            '342'
            '343'
            '344'
            '345'
            '346'
            '347'
            '348'
            '349'
            '350'
            '351'
            '352'
            '353'
            '354'
            '355'
            '356'
            '357'
            '358'
            '359'
            '360'
            '361'
            '362'
            '363'
            '364'
            '365'
            '366'
            '367'
            '368'
            '369'
            '370'
            '371'
            '372'
            '373'
            '374'
            '375'
            '376'
            '377'
            '378'
            '379'
            '380'
            '381'
            '382'
            '383'
            '384'
            '385'
            '386'
            '387'
            '388'
            '389'
            '390'
            '391'
            '392'
            '393'
            '394'
            '395'
            '396'
            '397'
            '398'
            '399'
            '400'
            '401'
            '402'
            '403'
            '404'
            '405'
            '406'
            '407'
            '408'
            '409'
            '410'
            '411'
            '412'
            '413'
            '414'
            '415'
            '416'
            '417'
            '418'
            '419'
            '420'
            '421'
            '422'
            '423'
            '424'
            '425'
            '426'
            '427'
            '428'
            '429'
            '430'
            '431'
            '432'
            '433'
            '434'
            '435'
            '436'
            '437'
            '438'
            '439'
            '440'
            '441'
            '442'
            '443'
            '444'
            '445'
            '446'
            '447'
            '448'
            '449'
            '450'
            '451'
            '452'
            '453'
            '454'
            '455'
            '456'
            '457'
            '458'
            '459'
            '460'
            '461'
            '462'
            '463'
            '464'
            '465'
            '466'
            '467'
            '468'
            '469'
            '470'
            '471'
            '472'
            '473'
            '474'
            '475'
            '476'
            '477'
            '478'
            '479'
            '480'
            '481'
            '482'
            '483'
            '484'
            '485'
            '486'
            '487'
            '488'
            '489'
            '490'
            '491'
            '492'
            '493'
            '494'
            '495'
            '496'
            '497'
            '498'
            '499'
            '500')
        end
      end
    end
  end
  object dlgColorcube: TColorDialog
    Left = 224
    Top = 528
  end
end
